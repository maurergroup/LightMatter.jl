function function_builder(sim::Simulation)
    laser=laser_factory(sim)
    if sim.athermalelectrons.EmbeddedAthEM == true
        sim_bg = build_backgroundTTM(sim)
        sys_embed = generate_expressions(sim,laser)
        sys_embed = Dict("$(key)_AthEM" => value for (key, value) in sys_embed)
        sys_bg = generate_expressions(sim_bg,laser)
        comb_sys = merge(sys_embed,sys_bg)
        return comb_sys
    else
        sys = generate_expressions(sim,laser)
        return sys
    end
end

function build_backgroundTTM(sim::Simulation)
    Tel = ElectronicTemperature(Enabled=true, Electron_PhononCoupling=true, Conductivity=sim.electronictemperature.Conductivity, 
                ElectronicHeatCapacity=sim.electronictemperature.ElectronicHeatCapacity, 
                ElectronPhononCouplingValue = sim.electronictemperature.ElectronPhononCouplingValue)
    Tph = PhononicTemperature(Enabled=true, Electron_PhononCoupling=true, Conductivity=sim.phononictemperature.Conductivity, 
                PhononicHeatCapacity=sim.electronictemperature.PhononicHeatCapacity)
    return build_Simulation(electronictemperature = Tel, phononictemperature=Tph)
end

function generate_expressions(sim::Simulation,laser::Expr)
    exprs = Dict{String,Union{Expr,Vector{Expr}}}()
    if sim.electronictemperature.Enabled == true
        merge!(exprs,Dict("Tel" => Lightmatter.electrontemperature_factory(sim,laser)))
    end
    if sim.phononictemperature.Enabled == true
        merge!(exprs,Dict("Tph" => Lightmatter.phonontemperature_factory(sim)))
    end
    if sim.athermalelectrons.Enabled == true
        merge!(exprs,Dict("fneq" => Lightmatter.athemdistribution_factory(sim,laser)))
        if sim.athermalelectrons.AthermalElectron_ElectronCoupling == true
            merge!(exprs,Dict("noe" => Lightmatter.athem_thermalelectronparticlechange(sim)))
            τee = electron_relaxationtime(sim::Simulation)
            merge!(exprs,Dict("relax" => :(Lightmatter.athem_electronelectronscattering(Tel,μ,sim,fneq,DOS,n,τee))))
        end
    end
    return exprs
end

function generate_initialconditions(sim::Simulation,initialtemps::Dict{String, <:Real})
    temp_u0 = Dict()
    if sim.athermalelectrons.Enabled == true
        merge!(temp_u0, Dict("fneq"=>zeros(sim.structure.dimension.length,length(sim.structure.egrid))))
    end
    if sim.athermalelectrons.AthermalElectron_ElectronCoupling == true 
        if typeof(sim.structure.DOS) == Vector{spl}
            no_part=zeros(sim.structure.dimension.length)
            for j in eachindex(sim.structure.dimension.grid)
                no_part[j] = get_thermalparticles(0.0,1e-32,sim.structure.DOS[j],sim.structure.egrid)
            end
        else
            no_part = fill(get_thermalparticles(0.0,1e-32,sim.structure.DOS,sim.structure.egrid),sim.structure.dimension.length)
        end
        merge!(temp_u0,Dict("noe" => no_part ))
    end
    if sim.electronictemperature.Enabled == true
        merge!(temp_u0, Dict("Tel" => fill(initialtemps["Tel"],sim.structure.dimension.length)))
    end
    if sim.phononictemperature.Enabled == true
        merge!(temp_u0, Dict("Tph" => fill(initialtemps["Tph"],sim.structure.dimension.length)))
    end
    namtup = NamedTuple((Symbol(key),value) for (key,value) in temp_u0)
    return NamedArrayPartition(namtup)
end

function generate_parameters(sim::Simulation,initialtemps::Dict{String, <:Real})
    p = (sim=sim,)
    if sim.electronictemperature.Conductivity == true
        p = (; p..., Tel_cond = zeros(sim.structure.dimension.length))
    end
    if sim.athermalelectrons.Conductivity == true
        p = (; p..., f_cond = zeros(sim.structure.dimension.length,length(sim.structure.egrid)))
    end
    if sim.phononictemperature.Conductivity == true
        p = (; p..., Tph_cond = zeros(sim.structure.dimension.length))
    end
    if sim.athermalelectrons.Conductivity == true && sim.athermalelectrons.AthermalElectron_ElectronCoupling == true
        p = (; p..., Δn = zeros(sim.structure.dimension.length))
    end
    if sim.athermalelectrons.Enabled == true && sim.athermalelectrons.AthermalElectron_ElectronCoupling == false
        if typeof(sim.structure.DOS) == Vector{spl}
            no_part=zeros(sim.structure.dimension.length)
            for j in eachindex(sim.structure.dimension.grid)
                no_part[j] = get_thermalparticles(0.0,1e-32,sim.structure.DOS[j],sim.structure.egrid)
            end
        else
            no_part = fill(get_thermalparticles(0.0,1e-32,sim.structure.DOS,sim.structure.egrid),sim.structure.dimension.length)
        end
        p = (; p..., Tel = initialtemps["Tel"], noe = no_part)
    end
    if sim.electronictemperature.Enabled == true && sim.electronictemperature.AthermalElectron_ElectronCoupling == false
        if typeof(sim.structure.DOS) == Vector{spl}
            no_part=zeros(sim.structure.dimension.length)
            for j in eachindex(sim.structure.dimension.grid)
                no_part[j] = get_thermalparticles(0.0,1e-32,sim.structure.DOS[j],sim.structure.egrid)
            end
        else
            no_part = fill(get_thermalparticles(0.0,1e-32,sim.structure.DOS,sim.structure.egrid),sim.structure.dimension.length)
        end
        p = (; p..., noe = no_part)
    end
    return p 
end

function simulation_construction(sys,sim::Simulation)
    if sim.structure.Elemental_System == 1
        return monometallic_system(sys,sim)
    else
        return antenna_reactor_system(sys,sim)
    end
end

function monometallic_system(sys,sim::Simulation)
    cond_exprs = conductivity_expressions(sim)
    loop_body = build_loopbody(sys,sim)
    expr_cond = Expr(:block,cond_exprs...)
    return quote
        println(t)
        $expr_cond
        Threads.@threads for i in 1:p.sim.structure.dimension.length
            $loop_body
        end
    end
end

function conductivity_expressions(sim::Simulation)
    cond_exprs = []
    if sim.electronictemperature.Conductivity == true
        push!(cond_exprs,:(Lightmatter.electrontemperature_conductivity!(u.Tel,p.sim.electronictemperature.κ,sim.structure.dimension.spacing,u.Tph,p.Tel_cond)))
    end
    if sim.phononictemperature.Conductivity == true
        push!(cond_exprs,:(Lightmatter.phonontemperature_conductivity!(u.Tph,p.sim.phononictemperature.κ,sim.structure.dimension.spacing,p.Tph_cond)))
    end
    if sim.athermalelectrons.Conductivity == true
        push!(cond_exprs,:(Lightmatter.electron_distribution_transport!(sim.athermalelectrons.v_g,u.fneq,p.f_cond,p.sim.structure.dimension.spacing)))
        if sim.athermalelectrons.AthermalElectron_ElectronCoupling == true
            push!(cond_exprs,:(n_cond = Lightmatter.thermal_particle_transport!(sim.athermalelectrons.v_g,p.mp.egrid,u.noe,p.dim)))
        end
    end 
    return cond_exprs
end


function build_loopbody(sys,sim::Simulation)
    exprs = Vector{Expr}(undef,0)
    push!(exprs,variable_renaming(sim))
    push!(exprs, :(μ = Lightmatter.find_chemicalpotential(n,Tel,DOS,sim.structure.egrid)))

    if sim.athermalelectrons.EmbeddedAthEM == true
        embedding = quote
            if i == 1
                relax_dis = $(sys["relax_AthEM"])
                du.noe = $(sys["noe_AthEM"])
                Δn = du.noe[1]
                du.fneq[1,:] .= $(sys["fneq_AthEM"])
                du.Tel[1] = $(sys["Tel_AthEM"])
                du.Tph[1] = $(sys["Tph_AthEM"])
            else
                du.Tel[i] = $(sys["Tel"])
                du.Tph[i] = $(sys["Tph"])
            end
        end
        push!(exprs,embedding)
    else
        if sim.electronictemperature.Enabled == true && sim.electronictemperature.AthermalElectron_ElectronCoupling == true
            push!(exprs,:(relax_dis = $(sys["relax"])))
            push!(exprs,:(du.noe[i,:] .= $(sys["noe"])))
            push!(exprs,:(Δn = du.noe[i]))
        end

        if sim.athermalelectrons.Enabled == true
            push!(exprs,:(du.fneq[i,:] .= $(sys["fneq"])))
        end
        if sim.electronictemperature.Enabled== true
            push!(exprs,:(du.Tel[i,:] .= $(sys["Tel"])))
        end
        if sim.phononictemperature.Enabled == true
            push!(exprs,:(du.Tph[i,:] .= $(sys["Tph"])))
        end
    end
    return Expr(:block, exprs...)
end

function variable_renaming(sim::Simulation)
    old_name = [:(p.sim)]
    new_name = [:sim]
    if typeof(sim.structure.DOS) == Vector{spl}
        push!(old_name, :(sim.structure.DOS[i]))
        push!(new_name, :DOS)
    else
        push!(old_name, :(sim.structure.DOS))
        push!(new_name, :DOS)
    end
    if sim.athermalelectrons.Enabled == true
        push!(old_name,:(u.fneq[i,:]))
        push!(new_name,:fneq)
        if sim.athermalelectrons.Conductivity == true
            push!(old_name,:(p.f_cond[i,:]))
            push!(new_name,:f_cond)
        end
        if sim.athermalelectrons.AthermalElectron_ElectronCoupling == false
            push!(old_name,:(p.Tel))
            push!(new_name,:Tel)
            push!(old_name, :(p.noe[i]))
            push!(new_name, :n)
        else 
            push!(old_name, :(u.noe[i]))
            push!(new_name, :n)
        end
    end
    if sim.electronictemperature.Enabled == true
        push!(old_name,:(u.Tel[i]))
        push!(new_name,:Tel)
        if sim.athermalelectrons.AthermalElectron_ElectronCoupling == false
            push!(old_name, :(p.noe[i]))
            push!(new_name, :n)
        end
        if sim.electronictemperature.Conductivity == true
            push!(old_name,:(p.Tel_cond[i]))
            push!(new_name,:Tel_cond)
        end
    end
    if sim.phononictemperature.Enabled == true
        push!(old_name,:(u.Tph[i]))
        push!(new_name,:Tph)
        if sim.phononictemperature.Conductivity == true
            push!(old_name,:(p.Tph_cond[i]))
            push!(new_name,:Tph_cond)
        end
    end
    old_name = Tuple(old_name)
    new_name = Tuple(new_name)
    #return Expr(:local, Expr(:(=), Expr(:tuple, new_name...), Expr(:tuple, old_name...)))
    #return Expr(:(=), Expr(:tuple,new_name...), Expr(:tuple,old_name...))
    assignments = [:(local $(lhs) = $(rhs)) for (lhs, rhs) in zip(new_name, old_name)]
    return quote
        $(assignments...)
    end
end



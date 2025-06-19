"""
    function_builder(sim::Simulation)
    
    Assembles the correct dictionary of equations for the subsystems that are propagated. Mainly used for
    seperating embedded and un-embedded methods

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Dictionary of subsystem names and their respective equation / expression
"""
function function_builder(sim::Simulation)
    laser=laser_factory(sim)
    if sim.athermalelectrons.EmbeddedAthEM == true
        sim_bg = build_backgroundTTM(sim)
        sys_embed = generate_expressions(sim, laser)
        sys_embed = Dict("$(key)_AthEM" => value for (key, value) in sys_embed)
        sys_bg = generate_expressions(sim_bg, laser)
        comb_sys = merge(sys_embed,sys_bg)
        return comb_sys
    else
        sys = generate_expressions(sim, laser)
        return sys
    end
end
"""
    build_backgroundTTM(sim::Simulation)
    
    Builds the structs for the TTM to be used in the sub-surface layers of the embedding AthEM nmethod

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Simulation struct setup to perform the background TTM for the embedded AthEM method
"""
function build_backgroundTTM(sim::Simulation)
    Tel = ElectronicTemperature(Enabled=true, Electron_PhononCoupling=true, Conductivity=sim.electronictemperature.Conductivity, 
                ElectronicHeatCapacity=sim.electronictemperature.ElectronicHeatCapacity, 
                ElectronPhononCouplingValue = sim.electronictemperature.ElectronPhononCouplingValue)
    Tph = PhononicTemperature(Enabled=true, Electron_PhononCoupling=true, Conductivity=sim.phononictemperature.Conductivity, 
                PhononicHeatCapacity=sim.electronictemperature.PhononicHeatCapacity)
    return build_Simulation(electronictemperature = Tel, phononictemperature=Tph)
end
"""
    generate_expressions(sim::Simulation, laser::Expr)
    
    Calculates and groups each subsystems expression in turn into a dictionary

    # Arguments
    - 'sim': Simulation settings and parameters
    - 'laser': Expression for the temporal evolution and spatial decay of the laser

    # Returns
    - Dictionary of subssytems and their respective expressions
"""
function generate_expressions(sim::Simulation, laser::Expr)
    exprs = Dict{String,Union{Expr,Vector{Expr}}}()
    if sim.electronictemperature.Enabled == true
        merge!(exprs,Dict("Tel" => Lightmatter.electrontemperature_factory(sim, laser)))
    end
    if sim.phononictemperature.Enabled == true
        merge!(exprs,Dict("Tph" => Lightmatter.phonontemperature_factory(sim)))
    end
    if sim.athermalelectrons.MagnetoTransport == true
        merge!(exprs,Dict("magneto" => Lightmatter.magnetotransport_equations(sim)))
    end
    if sim.athermalelectrons.Enabled == true
        merge!(exprs,Dict("fneq" => Lightmatter.athemdistribution_factory(sim, laser)))
        if sim.athermalelectrons.AthermalElectron_ElectronCoupling == true
            merge!(exprs,Dict("noe" => Lightmatter.athem_thermalelectronparticlechange(sim)))
            τee = electron_relaxationtime(sim::Simulation)
            merge!(exprs,Dict("relax" => :(Lightmatter.athem_electronelectronscattering(Tel, μ, sim, fneq, DOS, tot_n, $τee))))
        end
    end
    return exprs
end
"""
    generate_initialconditions(sim::Simulation, initialtemps::Dict{String, Float64})
    
    Generates the initial conditions (u0) NamedArrayPartition for the ODE 

    # Arguments
    - 'sim': Simulation settings and parameters
    - 'initialtemps': Dictionary containing initial temepratures for electronic and phononic baths

    # Returns
    - NamedArrayPartition containing the initial conditions of the simulation
"""
function generate_initialconditions(sim::Simulation, initialtemps::Dict{String, Float64})
    temp_u0 = Dict()
    if sim.athermalelectrons.Enabled == true
        merge!(temp_u0, Dict("fneq"=>zeros(sim.structure.dimension.length, length(sim.structure.egrid))))
    end
    if sim.athermalelectrons.AthermalElectron_ElectronCoupling == true 
        if typeof(sim.structure.DOS) == Vector{spl}
            if sim.structure.Elemental_System == 1
                no_part=zeros(sim.structure.dimension.length)
                for j in eachindex(sim.structure.dimension.grid)
                    no_part[j] = get_thermalparticles(0.0, 1e-32, sim.structure.DOS[j], sim.structure.egrid)
                end
            else
                no_part=zeros(sim.structure.dimension.length)
                for j in eachindex(sim.structure.dimension.grid)
                    x = mat_picker(sim.structure.dimension.grid[j], sim.structure.dimension.InterfaceHeight)
                    no_part[j] = get_thermalparticles(0.0, 1e-32, sim.structure.DOS[x], sim.structure.egrid)
                end
            end
        else
            no_part = fill(get_thermalparticles(0.0,1e-32, sim.structure.DOS, sim.structure.egrid), sim.structure.dimension.length)
        end
        merge!(temp_u0,Dict("noe" => no_part ))
    end
    if sim.electronictemperature.Enabled == true
        merge!(temp_u0, Dict("Tel" => fill(initialtemps["Tel"], sim.structure.dimension.length)))
    end
    if sim.phononictemperature.Enabled == true
        merge!(temp_u0, Dict("Tph" => fill(initialtemps["Tph"], sim.structure.dimension.length)))
    end
    namtup = NamedTuple((Symbol(key),value) for (key,value) in temp_u0)
    return NamedArrayPartition(namtup)
end
"""
    generate_initialconditions(sim::Simulation, initialtemps::Dict{String, Float64})
    
    Generates the parameters as a NamedTuple for the ODE 

    # Arguments
    - 'sim': Simulation settings and parameters
    - 'initialtemps': Dictionary containing initial temepratures for electronic and phononic baths

    # Returns
    - NamedTuple containing the parameters of the simulation
"""
function generate_parameters(sim::Simulation, initialtemps::Dict{String, Float64})
    if sim.structure.Elemental_System > 1
        p = (sim=sim,matsim=sim_seperation(sim))
    else
        p = (sim=sim,)
    end
    if sim.electronictemperature.Conductivity == true
        p = (; p..., Tel_cond = Vector{Float64}(undef,sim.structure.dimension.length))
    end
    if sim.athermalelectrons.Conductivity == true
        p = (; p..., f_cond = Matrix{Float64}(undef,sim.structure.dimension.length, length(sim.structure.egrid)))
    end
    if sim.phononictemperature.Conductivity == true
        p = (; p..., Tph_cond = Vector{Float64}(undef,sim.structure.dimension.length))
    end
    if sim.athermalelectrons.Conductivity == true && sim.athermalelectrons.AthermalElectron_ElectronCoupling == true
        p = (; p..., Δn = Vector{Float64}(undef,sim.structure.dimension.length))
    end
    if sim.athermalelectrons.MagnetoTransport == true
        p = (; p..., Δf_mt = Matrix{Float64}(undef,sim.structure.dimension.length, length(sim.structure.egrid)), g_k = Matrix{Float64}(undef,sim.structure.dimension.length, length(sim.structure.egrid)))
    end
    if sim.athermalelectrons.Enabled == true && sim.athermalelectrons.AthermalElectron_ElectronCoupling == false
        if typeof(sim.structure.DOS) == Vector{spl} && sim.structure.Elemental_System == 1
            no_part=zeros(sim.structure.dimension.length)
            for j in eachindex(sim.structure.dimension.grid)
                no_part[j] = get_thermalparticles(0.0, 1e-32, sim.structure.DOS[j], sim.structure.egrid)
            end
        elseif typeof(sim.structure.DOS) == Vector{spl} && sim.structure.Elemental_System > 1
            no_part=zeros(sim.structure.dimension.length)
            for j in eachindex(sim.structure.dimension.grid)
                mat = mat_picker(sim.structure.dimension.grid[j], sim.structure.dimension.InterfaceHeight)
                no_part[j] = get_thermalparticles(0.0, 1e-32, sim.structure.DOS[mat], sim.structure.egrid)
            end
        else
            no_part = fill(get_thermalparticles(0.0, 1e-32, sim.structure.DOS, sim.structure.egrid), sim.structure.dimension.length)
        end
        p = (; p..., Tel = initialtemps["Tel"], noe = no_part)
    end
    if sim.electronictemperature.Enabled == true && sim.electronictemperature.AthermalElectron_ElectronCoupling == false
        if typeof(sim.structure.DOS) == Vector{spl} && sim.structure.Elemental_System == 1
            no_part=zeros(sim.structure.dimension.length)
            for j in eachindex(sim.structure.dimension.grid)
                no_part[j] = get_thermalparticles(0.0, 1e-32, sim.structure.DOS[j], sim.structure.egrid)
            end
        elseif typeof(sim.structure.DOS) == Vector{spl} && sim.structure.Elemental_System > 1
            no_part=zeros(sim.structure.dimension.length)
            for j in eachindex(sim.structure.dimension.grid)
                mat = mat_picker(sim.structure.dimension.grid[j], sim.structure.dimension.InterfaceHeight)
                no_part[j] = get_thermalparticles(0.0, 1e-32, sim.structure.DOS[mat], sim.structure.egrid)
            end
        else
            no_part = fill(get_thermalparticles(0.0, 1e-32, sim.structure.DOS, sim.structure.egrid), sim.structure.dimension.length)
        end
        p = (; p..., noe = no_part)
    end
    return p 
end
"""
    simulation_construction(sys, sim::Simulation)
    
    Creates the expression block for the entire ODE function including multithreading.
    This function calls monometallic_system or antenna_reactor_system respectively.

    # Arguments
    - 'sys': Dictionary of expressions for each subssytem propagated
    - 'sim': Simulation settings and parameters

    # Returns
    - Quote block for the entire ODE problem
"""
function simulation_construction(sys, sim::Simulation)
    if sim.structure.Elemental_System == 1
        return monometallic_system(sys, sim)
    else
        return antenna_reactor_system(sys, sim)
    end
end
"""
    monometallic_system(sys, sim::Simulation)
    
    Creates the expression block for the entire ODE function including multithreading.
    This function is sepcfically for a monometallic system, Elemental_System == 1

    # Arguments
    - 'sys': Dictionary of expressions for each subssytem propagated
    - 'sim': Simulation settings and parameters

    # Returns
    - Quote block for the monometallic problem
"""
function monometallic_system(sys, sim::Simulation)
    expr_cond = conductivity_expressions(sim)
    loop_body = build_loopbody(sys, sim)
    return quote
        println(t)
        $expr_cond
        Threads.@threads for i in 1:p.sim.structure.dimension.length
            $loop_body
        end
    end
end
"""
    conductivity_expressions(sim::Simulation)
    
    Creates vector of expressions for the different conductivities occuring during the simulation

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Vector of expression for the conductivity of each subsytem if they are enabled
"""
function conductivity_expressions(sim::Simulation)
    cond_exprs = []
    if sim.electronictemperature.Conductivity == true
        push!(cond_exprs,:(Lightmatter.electrontemperature_conductivity!(u.Tel, p.sim.electronictemperature.κ, p.sim.structure.dimension.spacing, u.Tph, p.Tel_cond)))
    end
    if sim.phononictemperature.Conductivity == true
        push!(cond_exprs,:(Lightmatter.phonontemperature_conductivity!(u.Tph, p.sim.phononictemperature.κ, p.sim.structure.dimension.spacing, p.Tph_cond)))
    end
    if sim.athermalelectrons.Conductivity == true
        push!(cond_exprs,:(Lightmatter.electron_distribution_transport!(p.sim.athermalelectrons.v_g, u.fneq, p.f_cond, p.sim.structure.dimension.spacing)))
        if sim.athermalelectrons.AthermalElectron_ElectronCoupling == true
            push!(cond_exprs,:(n_cond = Lightmatter.thermal_particle_transport!(p.sim.athermalelectrons.v_g, p.sim.structure.egrid, u.noe, p.Δn, p.sim.structure.dimension.spacing)))
        end
    end 
    return Expr(:block,cond_exprs...)
end
"""
    build_loopbody(sys, sim::Simulation)
    
    Builds the multi-threaded section of the ODE problem

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Quote block for the multithreaded section of the ODE problem
"""
function build_loopbody(sys, sim::Simulation)
    exprs = Vector{Expr}(undef,0)
    push!(exprs,variable_renaming(sim))
    push!(exprs, :(μ = Lightmatter.find_chemicalpotential(n, Tel, DOS, sim.structure.egrid)))

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
            push!(exprs, :(tot_n = n + Lightmatter.get_noparticles(fneq, DOS, sim.structure.egrid)))
            push!(exprs,:(relax_dis = $(sys["relax"])))
            if sim.athermalelectrons.MagnetoTransport == true
                push!(exprs,:($(sys["magneto"])))
            end
            push!(exprs,:(du.noe[i] = $(sys["noe"])))
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
"""
    variable_renaming(sim::Simulation)
    
    Renames variables such as u.x or p.y to just x or y so that during expression evaluation they 
    are named correctly. 

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Expression for the variabble renaming to enter the top of the multithreaded loop
"""
function variable_renaming(sim::Simulation)
    old_name = [:(p.sim), :(p.sim.structure.tmp[i,:])]
    new_name = [:sim, :tmp]
    if typeof(sim.structure.DOS) == Vector{spl}
        push!(old_name, :(p.sim.structure.DOS[i]))
        push!(new_name, :DOS)
    else
        push!(old_name, :(p.sim.structure.DOS))
        push!(new_name, :DOS)
    end
    if sim.athermalelectrons.Enabled == true
        push!(old_name,:(u.fneq[i,:]))
        push!(new_name,:fneq)
        if sim.athermalelectrons.Conductivity == true
            push!(old_name,:(p.f_cond[i,:]))
            push!(new_name,:f_cond)
            if sim.athermalelectrons.AthermalElectron_ElectronCoupling == true
                push!(old_name,:(p.Δn[i]))
                push!(new_name,:n_cond)
            end
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
        if sim.athermalelectrons.MagnetoTransport == true
            push!(old_name, :(p.Δf_mt[i,:]))
            push!(old_name, :(p.g_k[i,:]))
            push!(new_name, :Δf_mt)
            push!(new_name, :g_k)
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
    assignments = [:(local $(lhs) = $(rhs)) for (lhs, rhs) in zip(new_name, old_name)]
    return quote
        $(assignments...)
    end
end
function function_builder(sim,las,dim)
    laser=laser_factory(las,dim)
    if sim.ParameterApprox.EmbeddingMethod == true && sim.Systems==SystemComponents(true,true,true)
        sim_bg = SimulationSettings(ParameterApproximation(true, true, true, true), Interaction(false, true), SystemComponents(true, true, false), true)
        sys_embed = generate_expressions(sim,laser,dim)
        sys_embed = Dict("$(key)_AthEM" => value for (key, value) in sys_embed)
        sys_bg = generate_expressions(sim_bg,laser,dim)
        comb_sys = merge(sys_embed,sys_bg)
        return comb_sys
    else
        sys = generate_expressions(sim,laser,dim)
        return sys
    end
end

function generate_expressions(sim,laser,dim)
    exprs = Dict{String,Union{Expr,Vector{Expr}}}()
    if sim.Systems.ElectronTemperature == true
        merge!(exprs,Dict("Tel" => electrontemperature_factory(sim,laser)))
    end
    if sim.Systems.PhononTemperature == true
        merge!(exprs,Dict("Tph" => phonontemperature_factory(sim)))
    end
    if sim.Systems.NonEqElectrons == true
        merge!(exprs,Dict("fneq" => athemdistribution_factory(sim,laser)))
        if sim.Interactions.ElectronElectron == true
            merge!(exprs,Dict("noe" => athem_electronparticlechange()))
            merge!(exprs,Dict("relax" => :(athem_electronelectronscattering(Tel,μ,mp,fneq,DOS,n))))
        end
    end
    return exprs
end

function generate_initialconditions(sim,mp,initialtemps,dim)
    if sim.Systems.NonEqElectrons == true
        if sim.Systems.ElectronTemperature == false
            temp_u0 = Dict("fneq"=>zeros(dim.length,length(mp.egrid)))
        elseif sim.Systems.PhononTemperature == false
            no_part=zeros(dim.length)
            for j in eachindex(dim.grid)
                no_part[j] = get_thermalparticles(0.0,1e-32,mp.DOS[j],8.617e-5,mp.egrid)
            end
            temp_u0 = Dict("noe"=>no_part, "fneq"=>zeros(dim.length,length(mp.egrid)), "Tel"=>fill(initialtemps["Tel"],dim.length))
        else
            no_part=zeros(dim.length)
            for j in eachindex(dim.grid)
                no_part[j] = get_thermalparticles(0.0,1e-32,mp.DOS[j],8.617e-5,mp.egrid)
            end
            temp_u0 = Dict("noe"=>no_part, "fneq"=>zeros(dim.length,length(mp.egrid)), "Tel"=>fill(initialtemps["Tel"],dim.length), "Tph"=>fill(initialtemps["Tph"],dim.length))
        end
    else
        temp_u0 = Dict("Tel"=>fill(initialtemps["Tel"],dim.length), "Tph"=>fill(initialtemps["Tel"],dim.length))
    end
    namtup = NamedTuple((Symbol(key),value) for (key,value) in temp_u0)
    return NamedArrayPartition(namtup)
end

function generate_parameters(sim,las,mp,initialtemps,dim)
    if sim.Systems.NonEqElectrons==true
        if sim.Systems.ElectronTemperature==false
            μ = find_chemicalpotential(mp.n0[1],initialtemps["Tel"],mp.DOS[1],cons.kB,mp.egrid)
            p = (las=las,mp=mp,dim=dim,Tel=initialtemps["Tel"],μ=μ)
        elseif sim.Systems.PhononTemperature == true
            p=(las=las,mp=mp,dim=dim,cond=zeros(dim.length))
        else
            p=(las=las,mp=mp,dim=dim)
        end
    else
        p=(las=las,mp=mp,dim=dim,noe=mp.n0,cond=zeros(dim.length))
    end
    return p
end

function simulation_construction(sys,sim)
    if sim.Systems.ElectronTemperature == true && sim.Systems.PhononTemperature == true
        expr_cond = :(Lightmatter.electrontemperature_conductivity!(u.Tel,p.dim,u.Tph,p.mp,p.cond))
    else
        expr_cond = :(nothing)
    end
    loop_body = build_loopbody(sys,sim)
    return quote
        println(t)
        $expr_cond
        Threads.@threads for i in 1:p.dim.length
            $loop_body
        end
    end
end

function build_loopbody(sys,sim)
    exprs = Vector{Expr}(undef,0)
    push!(exprs,variable_renaming(sim))
    if sim.Systems.ElectronTemperature == true && sim.Systems.NonEqElectrons == true
        push!(exprs,:(μ = Lightmatter.find_chemicalpotential(u.noe[i],u.Tel[i],mp.DOS[i],cons.kB,mp.egrid)))
    else 
        push!(exprs,:(μ = Lightmatter.find_chemicalpotential(mp.n0[i],Tel,mp.DOS[i],cons.kB,mp.egrid)))
    end

    if sim.ParameterApprox.EmbeddingMethod == true
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
        if sim.Systems.ElectronTemperature == true && sim.Systems.NonEqElectrons == true
            push!(exprs,:(relax_dis = $(sys["relax"])))
            push!(exprs,:(du.noe[i,:] .= $(sys["noe"])))
            push!(exprs,:(Δn = du.noe[i]))
        end

        if sim.Systems.NonEqElectrons == true
            push!(exprs,:(du.fneq[i,:] .= $(sys["fneq"])))
        end
        if sim.Systems.ElectronTemperature == true
            push!(exprs,:(du.Tel[i,:] .= $(sys["Tel"])))
        end
        if sim.Systems.PhononTemperature == true
            push!(exprs,:(du.Tph[i,:] .= $(sys["Tph"])))
        end
    end
    return Expr(:block, exprs...)
end

function variable_renaming(sim)
    old_name=[:(p.mp.DOS[i]),:(p.mp),:(p.las),:(p.dim)]
    new_name=[:DOS,:mp,:las,:dim]
    if sim.Systems.NonEqElectrons == true
        push!(old_name,:(u.fneq[i,:]))
        push!(new_name,:fneq)
        if sim.Systems.ElectronTemperature == false
            push!(old_name,:(p.Tel))
            push!(new_name,:Tel)
        else 
            push!(old_name, :(u.noe[i]))
            push!(new_name, :n)
        end
    end
    if sim.Systems.ElectronTemperature == true
        push!(old_name,:(u.Tel[i]))
        push!(new_name,:Tel)
        if sim.Systems.PhononTemperature == true
            push!(old_name,:(p.cond))
            push!(new_name,:cond)
        end
    end
    if sim.Systems.PhononTemperature == true
        push!(old_name,:(u.Tph[i]))
        push!(new_name,:Tph)
        push!(old_name,:(p.cond[i]))
        push!(new_name,:cond)
    end
    old_name = Tuple(old_name)
    new_name = Tuple(new_name)
    return Expr(:(=), Expr(:tuple,new_name...), Expr(:tuple,old_name...))
end
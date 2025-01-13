function function_builder(sim,las,dim;sys_output=false)
    laser=laser_factory(las,dim)
    if sim.ParameterApprox.EmbeddingMethod == true && sim.Systems==SystemComponents(true,true,true)
        sim_bg = SimulationSettings(ParameterApproximation(true, true, true, true), Interaction(false, true), SystemComponents(true, true, false), true)
        sys_embed,key_list_embed = generate_expressions(sim,laser,dim)
        sys_embed = Dict("$(key)_AthEM" => value for (key, value) in sys_embed)
        key_list_embed = ["$(s)_AthEM" for s in key_list_embed]
        sys_bg,key_list_bg = generate_expressions(sim_bg,laser,dim)
        args_embed = generate_arguments(sim)
        args_embed = Dict("$(key)_AthEM" => value for (key, value) in args_embed)
        args_bg = generate_arguments(sim_bg)

        scalar_functions(sys_embed,key_list_embed,args_embed)
        multithread_functions(sys_bg,key_list_bg,args_bg)
        if sys_output==true
            return merge(sys_embed,sys_bg),vcat(key_list_embed,key_list_bg)
        else
            return vcat(key_list_embed,key_list_bg)
        end
    else
        sys,key_list = generate_expressions(sim,laser,dim)
        args = generate_arguments(sim)
        if typeof(dim) == Homogeneous
            scalar_functions(sys,key_list,args)
        else
            multithread_functions(sys,key_list,args)
        end
        if sys_output==true
            return sys,key_list
        else
            return key_list
        end
    end
end

function generate_expressions(sim,laser,dim)
    exprs = Dict{String,Union{Expr,Vector{Expr}}}()
    if sim.Systems.ElectronTemperature == true
        merge!(exprs,Dict("Tel" => electrontemperature_factory(sim,laser,dim)))
    end
    if sim.Systems.PhononTemperature == true
        merge!(exprs,Dict("Tph" => phonontemperature_factory(sim)))
    end
    if sim.Systems.NonEqElectrons == true
        merge!(exprs,Dict("fneq" => athemdistribution_factory(sim,laser)))
        if sim.Interactions.ElectronElectron == true
            merge!(exprs,Dict("noe" => athem_electronparticlechange()))
            merge!(exprs,Dict("relax" => athem_electronelectronscattering()))
        end
    end
    return exprs,collect(keys(exprs))    
end

function generate_arguments(sim::SimulationSettings)
    args = Dict{String,Tuple}()
    if sim.Systems.ElectronTemperature == true
        if sim.Interactions.ElectronElectron == true
            if sim.Interactions.ElectronPhonon == true
                merge!(args,Dict("Tel" => ((:Tel,Float64),(:Tph,Float64),(:mp,MaterialParameters),(:DOS,spl),(:cons,Constants),(:μ,Float64),(:relax_dis,Vector{Float64}),(:Δn,Float64),(:cond,Float64))))
            else
                merge!(args,Dict("Tel" => ((:Tel,Float64),(:mp,MaterialParameters),(:DOS,spl),(:cons,Constants),(:μ,Float64),(:relax_dis,Vector{Float64}),(:Δn,Float64),(:cond,Float64))))
            end
        elseif sim.Interactions.ElectronPhonon == true
            merge!(args,Dict("Tel" => ((:Tel,Float64),(:Tph,Float64),(:mp,MaterialParameters),(:DOS,spl),(:cons,Constants),(:las,Laser),(:μ,Float64),(:t,Float64),(:cond,Float64),(:dim,Dimension))))
        else
            merge!(args,Dict("Tel" => ((:Tel,Float64),(:mp,MaterialParameters),(:DOS,spl),(:cons,Constants),(:las,Laser),(:μ,Float64),(:t,Float64),(:cond,Float64))))
        end
    end


    if sim.Systems.PhononTemperature == true
        if sim.Systems.NonEqElectrons == true
            merge!(args,Dict("Tph" => ((:Tel,Float64),(:Tph,Float64),(:mp,MaterialParameters),(:DOS,spl),(:cons,Constants),(:μ,Float64),(:fneq,Vector{Float64}))))
        else
            merge!(args,Dict("Tph" => ((:Tel,Float64),(:Tph,Float64),(:mp,MaterialParameters),(:DOS,spl),(:cons,Constants),(:μ,Float64))))
        end
    end

    if sim.Systems.NonEqElectrons == true
        if sim.Interactions.ElectronElectron == true
            merge!(args,Dict("fneq" => ((:fneq,Vector{Float64}),(:Tel,Float64),(:mp,MaterialParameters),(:DOS,spl),(:cons,Constants),(:las,Laser)
                                ,(:μ,Float64),(:t,Float64),(:dim,Dimension),(:relax_dis,Vector{Float64}))))
            merge!(args,Dict("noe" => ((:relax_dis,Vector{Float64}),(:μ,Float64),(:mp,MaterialParameters),(:DOS,spl))))
            merge!(args,Dict("relax" => ((:Tel,Float64),(:fneq,Vector{Float64}),(:n,Float64),(:μ,Float64),(:mp,MaterialParameters),(:DOS,spl),(:cons,Constants))))
        else
            merge!(args,Dict("fneq" => ((:fneq,Vector{Float64}),(:Tel,Float64),(:mp,MaterialParameters),(:DOS,spl),(:cons,Constants),(:las,Laser)
                                        ,(:μ,Float64),(:t,Float64),(:dim,Dimension))))
        end
    end

    return args
end

function generate_initialconditions(key_list,mp,initialtemps,dim)
    temp_u = ()
    for i in key_list
        if i =="Tel"
            temp_u = (temp_u...,fill(initialtemps[i],dim.length))
        elseif i == "Tph"
            temp_u = (temp_u...,fill(initialtemps[i],dim.length))
        elseif i == "fneq"
            temp_u = (temp_u...,zeros(dim.length,length(mp.egrid)))
        elseif i == "noe"
            no_part=zeros(dim.length)
            for j in eachindex(dim.grid)
                no_part[j] = get_thermalparticles(0.0,1e-32,mp.DOS[j],8.617e-5,mp.egrid)
            end
            temp_u = (temp_u...,no_part)
        elseif i == "fneq_AthEM"
            temp_u = (temp_u...,zeros(1,length(mp.egrid)))
        elseif i == "noe_AthEM"
            no_part=zeros(dim.length)
            for j in eachindex(dim.grid)
                no_part[j] = get_thermalparticles(0.0,1e-32,mp.DOS[j],8.617e-5,mp.egrid)
            end
            temp_u = (temp_u...,no_part)
        end
    end
    return ArrayPartition(temp_u)
end

function generate_parameters(sim,mp,cons,las,initialtemps,dim)
    if sim.Systems.NonEqElectrons==true
        if sim.Systems.ElectronTemperature==false
            μ = find_chemicalpotential(mp.n0,initialtemps["Tel"],mp.DOS[1],cons.kB,mp.egrid)
            return [las,mp,cons,dim,initialtemps["Tel"],μ]
        else
            if sim.ParameterApprox.EmbeddingMethod == true
                return [las,mp,cons,dim,zeros(dim.length),zeros(dim.length),zeros(1,length(mp.egrid))]
            else
                return return [las,mp,cons,dim,zeros(dim.length),zeros(dim.length),zeros(dim.length,length(mp.egrid))]
            end
        end
    else
        return [las,mp,cons,dim,mp.n0,zeros(dim.length),zeros(dim.length)]
    end
end

function scalar_functions(sys,key_list,args) #Generates Scalar functions
    for i in key_list #Cycles through the different functions
        vars = collect(x for (x,y) in args[i])
        if i == "fneq_AthEM"
            push!(vars,:z)
        end
        make_function(vars,sys[i],Symbol(i*"_scal"))
        new_args = []
        scalar_args = []
        for j in args[i]
            if j[2] == Float64 && j[1] != :t
                push!(new_args,:($(j[1])[1]))
                push!(scalar_args,j[1])
            elseif j[2] == Vector{Float64} 
                push!(new_args,:($(j[1])[:]))
                push!(scalar_args,j[1])
            elseif j[2] == spl
                push!(new_args,:($(j[1])[1]))
                push!(scalar_args,j[1])
            else
                push!(new_args,j[1])
                push!(scalar_args,j[1])
            end
        end
        if i == "fneq_AthEM"
            push!(new_args,1)
        end
        expr = scalar_expr(Symbol(i*"_scal"),:du,new_args)
        push!(scalar_args,:du)
        make_function(scalar_args,expr,Symbol(i*"_func"))
    end
end

function multithread_functions(sys,key_list,args)
    for i in key_list
        scal_args = collect(x for (x,y) in args[i])
        scal_args = push!(scal_args,:z)
        make_function(scal_args,sys[i],Symbol(i*"_scal"))
        new_args = []
        parallel_args = []
        for j in args[i]
            if j[2] == Float64 && j[1] != :t
                push!(new_args,:($(j[1])[i]))
                push!(parallel_args,j[1])
            elseif j[2] == Vector{Float64} 
                push!(new_args,:($(j[1])[i,:]))
                push!(parallel_args,j[1])
            elseif j[2] == spl
                push!(new_args,:($(j[1])[i]))
                push!(parallel_args,j[1])
            else
                push!(new_args,j[1])
                push!(parallel_args,j[1])
            end
        end
        expr = multithreaded_expr(Symbol(i*"_scal"),:du,new_args)
        push!(parallel_args,:du)
        make_function(parallel_args,expr,Symbol(i*"_func"))
    end
end  

function make_function(args, expr::Expr,name::Symbol)
    function_expr = Expr(:function, Expr(:call, name, args...), expr)
    eval(function_expr)
end

function multithreaded_expr(func, result::Symbol, vars)
    return quote
        Threads.@threads for i in eachindex(du[:,1])
            @inbounds $result[i,:] .= $func($(vars...),i)
        end
    end
end

function scalar_expr(func,result,vars)
    return quote
        $result[:] .= $func($(vars...))
    end
end
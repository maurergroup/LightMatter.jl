function function_builder()
    sim,mp,las,dim,cons = setup()
    laser=laser_factory(las,dim)
    sys,key_list = generate_expressions(sim,laser,dim)
    args = generate_arguments(sim)
    if typeof(dim) == Homogeneous
        scalar_functions(sys,key_list,args)
    else
        multithread_functions(sys,key_list,args)
    end
    return key_list
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
        merge!(exprs,Dict("fneq" => athemdistribution_factory(sim,laser,dim)))
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
                merge!(args,Dict("Tel" => ((:Tel,Real),(:Tph,Real),(:mp,MaterialParameters),(:cons,Constants),(:μ,Real),(:relax_dis,Vector{<:Real}),(:Δn,Real),(:cond,Real))))
            else
                merge!(args,Dict("Tel" => ((:Tel,Real),(:mp,MaterialParameters),(:cons,Constants),(:μ,Real),(:relax_dis,Vector{<:Real}),(:Δn,Real),(:cond,Real))))
            end
        elseif sim.Interactions.ElectronPhonon == true
            merge!(args,Dict("Tel" => ((:Tel,Real),(:Tph,Real),(:mp,MaterialParameters),(:cons,Constants),(:las,Laser),(:μ,Real),(:t,Real),(:cond,Real),(:dim,Dimension))))
        else
            merge!(args,Dict("Tel" => ((:Tel,Real),(:mp,MaterialParameters),(:cons,Constants),(:las,Laser),(:μ,Real),(:t,Real),(:cond,Real))))
        end
    end


    if sim.Systems.PhononTemperature == true
        if sim.Systems.NonEqElectrons == true
            merge!(args,Dict("Tph" => ((:Tel,Real),(:Tph,Real),(:mp,MaterialParameters),(:cons,Constants),(:μ,Real),(:relax_Tphdis,Vector{<:Real}))))
        else
            merge!(args,Dict("Tph" => ((:Tel,Real),(:Tph,Real),(:mp,MaterialParameters),(:cons,Constants),(:μ,Real))))
        end
    end

    if sim.Systems.NonEqElectrons == true
        if sim.Interactions.ElectronElectron == true
            merge!(args,Dict("fneq" => ((:fneq,Vector{<:Real}),(:Tel,Real),(:mp,MaterialParameters),(:cons,Constants),(:las,Laser)
                                        ,(:μ,Real),(:t,Real),(:dim,Dimension),(:relax_dis,Vector{<:Real}))))
            merge!(args,Dict("noe" => ((:relax_dis,Vector{<:Real}),(:μ,Real),(:mp,MaterialParameters))))
            merge!(args,Dict("relax" => ((:Tel,Real),(:fneq,Vector{<:Real}),(:n,Real),(:μ,Real),(:mp,MaterialParameters),(:cons,Constants))))
        else
            merge!(args,Dict("fneq" => ((:fneq,Vector{<:Real}),(:Tel,Real),(:mp,MaterialParameters),(:cons,Constants),(:las,Laser)
                                        ,(:μ,Real),(:t,Real),(:dim,Dimension))))
        end
    end

    return args
end

function generate_initalconditions(key_list,mp,cons,initialtemps,dim)
    temp_u = ()
    for i in key_list
        if i =="Tel"
            temp_u = (temp_u...,fill(initialtemps[i],dim.length))
        elseif i == "Tph"
            temp_u = (temp_u...,fill(initialtemps[i],dim.length))
        elseif i == "fneq"
            temp_u = (temp_u...,zeros(dim.length,length(mp.egrid)))
        elseif i == "noe"
            temp_u = (temp_u...,fill(get_thermalparticles(mp.μ,1e-16,mp.DOS,cons.kB,mp.FE),dim.length))
        end
    end
    return ArrayPartition(temp_u)
end

function generate_parameters(sim,mp,cons,las,initialtemps,dim)
    if sim.Systems.NonEqElectrons==true
        if sim.Systems.ElectronTemperature==false
            n = get_thermalparticles(mp.μ,1e-16,mp.DOS,cons.kB,mp.FE)
            μ = find_chemicalpotential(n,initialtemps["Tel"],mp.DOS,cons.kB,mp.FE)
            return (las,mp,cons,dim,initaltemps["Tel"],n,μ)
        else
            return (las,mp,cons,dim)
        end
    else
        n = get_thermalparticles(0.0,1e-16,mp.DOS,cons.kB,mp.FE)
        return (las,mp,cons,dim,n)
    end
end

function scalar_functions(sys,key_list,args)
    for i in key_list
        make_function(args[i],sys[i],Symbol(i*"_scal"))
        new_args = ()
        scalar_args = ()
        for j in args[i]
            if j[2] == Real && j[1] != :t
                new_args=(new_args...,:($(j[1])[1]))
                scalar_args=(scalar_args...,(j[1],Vector{<:Real}))
            elseif j[2] == Vector{<:Real} 
                new_args=(new_args...,:($(j[1])[1,:]))
                scalar_args=(scalar_args...,(j[1],Matrix{<:Real}))
            else
                new_args=(new_args...,j[1])
                scalar_args=(scalar_args...,j)
            end
        end
        expr = scalar_expr(Symbol(i*"_scal"),:du,new_args)
        if i == "relax" || i == "fneq"
            scalar_args = (scalar_args...,(:du,Matrix{<:Real}))
        else
            scalar_args = (scalar_args...,(:du,Vector{<:Real}))
        end
        make_function(scalar_args,expr,Symbol(i))
    end
end

function multithread_functions(sys,key_list,args)
    for i in key_list
        scal_args = (args[i]...,(:i,Int))
        make_function(scal_args,sys[i],Symbol(i*"_scal"))
        new_args = ()
        parallel_args = ()
        for j in args[i]
            if j[2] == Real && j[1] != :t
                new_args=(new_args...,:($(j[1])[i]))
                parallel_args=(parallel_args...,(j[1],Vector{<:Real}))
            elseif j[2] == Vector{<:Real} 
                new_args=(new_args...,:($(j[1])[i,:]))
                parallel_args=(parallel_args...,(j[1],Matrix{<:Real}))
            else
                new_args=(new_args...,j[1])
                parallel_args=(parallel_args...,j)
            end
        end
        expr = multithreaded_expr(Symbol(i*"_scal"),:du,new_args)
        if i == "relax" || i == "fneq"
            parallel_args = (parallel_args...,(:du,Matrix{<:Real}))
        else
            parallel_args = (parallel_args...,(:du,Vector{<:Real}))
        end
        make_function(parallel_args,expr,Symbol(i))
    end
end  

function make_function(typed_vars::Tuple, expr::Expr,name::Symbol)
    args = [Expr(:(::), var, typ) for (var, typ) in typed_vars]
    function_expr = Expr(:function, Expr(:call, name, args...), expr)
    eval(function_expr)
    eval(function_expr.args[1].args[1])
end

function multithreaded_expr(func, result::Symbol, vars)
    return quote
        Threads.@threads for i in 1:length($result[:,1])
            @inbounds $result[i] = $func($(vars...),i)
        end
    end
end

function scalar_test(func, result::Symbol, vars)
    return quote
        for i in 1:length($result[:,1])
            @inbounds $result[i] = $func($(vars...),i)
        end
    end
end

function scalar_expr(func,result,vars)
    return quote
        $result[:] .= $func($(vars...))
    end
end
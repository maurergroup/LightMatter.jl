function function_builder(sim,las,dim)
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
                merge!(args,Dict("Tel" => ((:Tel,Float64),(:Tph,Float64),(:mp,MaterialParameters),(:cons,Constants),(:μ,Float64),(:relax_dis,Vector{Float64}),(:Δn,Float64),(:cond,Float64))))
            else
                merge!(args,Dict("Tel" => ((:Tel,Float64),(:mp,MaterialParameters),(:cons,Constants),(:μ,Float64),(:relax_dis,Vector{Float64}),(:Δn,Float64),(:cond,Float64))))
            end
        elseif sim.Interactions.ElectronPhonon == true
            merge!(args,Dict("Tel" => ((:Tel,Float64),(:Tph,Float64),(:mp,MaterialParameters),(:cons,Constants),(:las,Laser),(:μ,Float64),(:t,Float64),(:cond,Float64),(:dim,Dimension))))
        else
            merge!(args,Dict("Tel" => ((:Tel,Float64),(:mp,MaterialParameters),(:cons,Constants),(:las,Laser),(:μ,Float64),(:t,Float64),(:cond,Float64))))
        end
    end


    if sim.Systems.PhononTemperature == true
        if sim.Systems.NonEqElectrons == true
            merge!(args,Dict("Tph" => ((:Tel,Float64),(:Tph,Float64),(:mp,MaterialParameters),(:cons,Constants),(:μ,Float64),(:fneq,Vector{Float64}))))
        else
            merge!(args,Dict("Tph" => ((:Tel,Float64),(:Tph,Float64),(:mp,MaterialParameters),(:cons,Constants),(:μ,Float64))))
        end
    end

    if sim.Systems.NonEqElectrons == true
        if sim.Interactions.ElectronElectron == true
            merge!(args,Dict("fneq" => ((:fneq,Vector{Float64}),(:Tel,Float64),(:mp,MaterialParameters),(:cons,Constants),(:las,Laser)
                                ,(:μ,Float64),(:t,Float64),(:dim,Dimension),(:relax_dis,Vector{Float64}))))
            merge!(args,Dict("noe" => ((:relax_dis,Vector{Float64}),(:μ,Float64),(:mp,MaterialParameters))))
            merge!(args,Dict("relax" => ((:Tel,Float64),(:fneq,Vector{Float64}),(:n,Float64),(:μ,Float64),(:mp,MaterialParameters),(:cons,Constants))))
        else
            merge!(args,Dict("fneq" => ((:fneq,Vector{Float64}),(:Tel,Float64),(:mp,MaterialParameters),(:cons,Constants),(:las,Laser)
                                        ,(:μ,Float64),(:t,Float64),(:dim,Dimension))))
        end
    end

    return args
end

function generate_initalconditions(key_list,mp,initialtemps,dim)
    temp_u = ()
    for i in key_list
        if i =="Tel"
            temp_u = (temp_u...,fill(initialtemps[i],dim.length))
        elseif i == "Tph"
            temp_u = (temp_u...,fill(initialtemps[i],dim.length))
        elseif i == "fneq"
            temp_u = (temp_u...,zeros(dim.length,length(mp.egrid)))
        elseif i == "noe"
            temp_u = (temp_u...,fill(mp.n0,dim.length))
        end
    end
    return ArrayPartition(temp_u)
end

function generate_parameters(sim,mp,cons,las,initialtemps,dim)
    if sim.Systems.NonEqElectrons==true
        if sim.Systems.ElectronTemperature==false
            μ = find_chemicalpotential(mp.n0,initialtemps["Tel"],mp.DOS,cons.kB,mp.FE,mp.n0)
            return [las,mp,cons,dim,initialtemps["Tel"],μ]
        else
            return [las,mp,cons,dim,zeros(dim.length),zeros(dim.length),zeros(dim.length,length(mp.egrid))]
        end
    else
        return [las,mp,cons,dim,mp.n0,zeros(dim.length),zeros(dim.length)]
    end
end

function scalar_functions(sys,key_list,args)
    for i in key_list
        make_function(args[i],sys[i],Symbol(i*"_scal"))
        new_args = ()
        scalar_args = ()
        for j in args[i]
            if j[2] == Float64 && j[1] != :t
                new_args=(new_args...,:($(j[1])[1]))
                scalar_args=(scalar_args...,(j[1],Vector{Float64}))
            elseif j[2] == Vector{Float64} 
                new_args=(new_args...,:($(j[1])[:]))
                scalar_args=(scalar_args...,(j[1],Matrix{Float64}))
            else
                new_args=(new_args...,j[1])
                scalar_args=(scalar_args...,j)
            end
        end
        expr = scalar_expr(Symbol(i*"_scal"),:du,new_args)
        if i == "relax" || i == "fneq"
            scalar_args = (scalar_args...,(:du,Matrix{Float64}))
        else
            scalar_args = (scalar_args...,(:du,Vector{Float64}))
        end
        make_function(scalar_args,expr,Symbol(i*"_func"))
    end
end

function multithread_functions(sys,key_list,args)
    for i in key_list
        scal_args = (args[i]...,(:z,Int))
        make_function(scal_args,sys[i],Symbol(i*"_scal"))
        new_args = ()
        parallel_args = ()
        for j in args[i]
            if j[2] == Float64 && j[1] != :t
                new_args=(new_args...,:($(j[1])[i]))
                parallel_args=(parallel_args...,(j[1],Vector{Float64}))
            elseif j[2] == Vector{Float64} 
                new_args=(new_args...,:($(j[1])[i,:]))
                parallel_args=(parallel_args...,(j[1],Matrix{Float64}))
            else
                new_args=(new_args...,j[1])
                parallel_args=(parallel_args...,j)
            end
        end
        expr = multithreaded_expr(Symbol(i*"_scal"),:du,new_args)
        if i == "relax" || i == "fneq"
            parallel_args = (parallel_args...,(:du,Matrix{Float64}))
        else
            parallel_args = (parallel_args...,(:du,Vector{Float64}))
        end
        make_function(parallel_args,expr,Symbol(i*"_func"))
    end
end  

function make_function(typed_vars::Tuple, expr::Expr,name::Symbol)
    args = [Expr(:(::), var, typ) for (var, typ) in typed_vars]
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
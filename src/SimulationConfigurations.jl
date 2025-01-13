function run_simulation(key_list,initialtemps,tspan,sim,mp,las,dim;save=2.0,tolerance=1e-4,max_step=0.1)
    u0 = generate_initialconditions(key_list,mp,initialtemps,dim)
    p = generate_parameters(sim,mp,cons,las,initialtemps,dim)
    return run_dynamics(key_list,p,u0,tspan,save,tolerance,max_step)
end

function run_dynamics(key_list,p,u0,tspan,save,tolerance,max_step)
    if length(key_list) == 1
        return run_dynamics1(p,u0,tspan,save,tolerance,max_step)
    elseif length(key_list) == 2
        return run_dynamics2(p,u0,tspan,save,tolerance,max_step)
    elseif length(key_list) == 3
        return sol = run_dynamics3(p,u0,tspan,save,tolerance,max_step)
    elseif length(key_list) == 5
        return sol = run_dynamics4(p,u0,tspan,save,tolerance,max_step)
    elseif length(key_list) == 7
        return sol = run_dynamics7(p,u0,tspan,save,tolerance,max_step)
    end
end

function run_dynamics1(p,u0,tspan,save,tolerance,max_step)
    println("Running small time span to precompile dynamics")
    solve(ODEProblem(EHP_simulation,u0,(0.0,0.1),p),Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmax=max_step)
    println("Running main dynamics")
    prob=ODEProblem(EHP_simulation,u0,tspan,p)
    sol = solve(prob,Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmax=max_step)
    return sol
end

function run_dynamics2(p,u0,tspan,save,tolerance,max_step)
    println("Running small time span to precompile dynamics")
    solve(ODEProblem(TTM_simulation,u0,(0.0,0.1),p),Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmax=max_step)
    println("Running main dynamics")
    prob=ODEProblem(TTM_simulation,u0,tspan,p)
    sol = solve(prob,Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmax=max_step)
    return sol
end

function run_dynamics3(p,u0,tspan,save,tolerance,max_step)
    println("Running small time span to precompile dynamics")
    solve(ODEProblem(eeAthEM_simulation,u0,(0.0,0.1),p),Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmax=max_step)
    println("Running main dynamics")
    prob=ODEProblem(eeAthEM_simulation,u0,tspan,p)
    sol = solve(prob,Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmax=max_step)
    return sol
end

function run_dynamics4(p,u0,tspan,save,tolerance,max_step)
    println("Running small time span to precompile dynamics")
    solve(ODEProblem(FullAthEM_simulation,u0,(0.0,0.1),p),Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmax=max_step)
    println("Running main dynamics")
    prob=ODEProblem(FullAthEM_simulation,u0,tspan,p)
    sol = solve(prob,Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmax=max_step)
    return sol
end

function run_dynamics7(p,u0,tspan,save,tolerance,max_step)
    println("Running small time span to precompile dynamics")
    solve(ODEProblem(embedded_AthEM_simulation,u0,(0.0,0.1),p),Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmax=max_step)
    println("Running main dynamics")
    prob=ODEProblem(embedded_AthEM_simulation,u0,tspan,p)
    sol = solve(prob,Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmax=max_step)
    return sol
end

function TTM_simulation(du,u,p,t)
    println(t)
    @assert all(!isnan(x) for x in u)
    electrontemperature_conductivity(u.x[1],p[4],u.x[2],p[2],p[7])
    Threads.@threads for i in eachindex(u.x[1])
        @inbounds p[6][i] = find_chemicalpotential(p[5],u.x[1][i],p[2].DOS[i],p[3].kB,p[2].egrid)
    end
    Tel_func(u.x[1],u.x[2],p[2],p[2].DOS,p[3],p[1],p[6],t,p[7],p[4],du.x[1])
    Tph_func(u.x[1],u.x[2],p[2],p[2].DOS,p[3],p[6],du.x[2])
    nothing
end

function eeAthEM_simulation(du,u,p,t)
    println(t) 
    #n = 1, Tel = 2, fneq = 3 

    μ = find_chemicalpotential.(u.x[1],u.x[2],Ref(p[2].DOS[1]),p[3].kB,p[2].egrid)
    relax_dis = zeros(length(u.x[3]),1)
    relax_func(u.x[2],u.x[3],u.x[1],μ,p[2],p[2].DOS,p[3],relax_dis)
    noe_func(relax_dis,μ,p[2],p[2].DOS,du.x[1])
    Tel_func(u.x[2],p[2],p[2].DOS,p[3],μ,relax_dis,du.x[1],[0.0],du.x[2])
    fneq_func(u.x[3],u.x[2],p[2],p[2].DOS,p[3],p[1],μ,t,p[4],relax_dis,du.x[3])
    nothing
end

function FullAthEM_simulation(du,u,p,t)
    println(t)

    electrontemperature_conductivity(u.x[1],p[4],u.x[3],p[2],p[6])
    Threads.@threads for i in eachindex(p[5])
        p[5][i] = find_chemicalpotential(u.x[2][i],u.x[1][i],p[2].DOS[i],p[3].kB,p[2].egrid)
    end
    relax_func(u.x[1],u.x[4],u.x[2],p[5],p[2],p[2].DOS,p[3],p[7])
    noe_func(p[7],p[5],p[2],p[2].DOS,du.x[2])
    Tel_func(u.x[1],u.x[3],p[2],p[2].DOS,p[3],p[5],p[7],du.x[2],p[6],du.x[1])
    Tph_func(u.x[1],u.x[3],p[2],p[2].DOS,p[3],p[5],u.x[4],du.x[3])
    fneq_func(u.x[4],u.x[1],p[2],p[2].DOS,p[3],p[1],p[5],t,p[4],p[7],du.x[4])
    nothing
end

function EHP_simulation(du,u,p,t)
    println(t)
    fneq_func(u.x[1],[p[5]],p[2],p[2].DOS,p[3],p[1],[p[6]],t,p[4],du.x[1])
    nothing
end

function embedded_AthEM_simulation(du,u,p,t)
    println(t)
    #electrontemperature_conductivity(u.x[3],p[4],u.x[4],p[2],p[6])
    electrontemperature_conductivity(u.x[3][2:end],p[4],u.x[4][2:end],p[2],view(p[6],2:p[4].length))
    p[6][1] = embedded_AthEM_conductivity(u.x[3],u.x[4],p[4],p[2])
    Threads.@threads for i in eachindex(p[5])
        p[5][i] = find_chemicalpotential(u.x[1][i],u.x[3][i],p[2].DOS[i],p[3].kB,p[2].egrid)
    end

    relax_AthEM_func(u.x[3][1],u.x[2][1,:],u.x[1][1],p[5][1],p[2],p[2].DOS,p[3],p[7])
    noe_AthEM_func(p[7],p[5][1],p[2],p[2].DOS,view(du.x[1],1))
    Tel_AthEM_func(u.x[3][1],u.x[4][1],p[2],p[2].DOS,p[3],p[5][1],p[7],du.x[1][1],p[6][1],view(du.x[3],1))
    Tph_AthEM_func(u.x[3][1],u.x[4][1],p[2],p[2].DOS,p[3],p[5][1],u.x[2][1,:],view(du.x[4],1))
    fneq_AthEM_func(u.x[2][1,:],u.x[3][1],p[2],p[2].DOS,p[3],p[1],p[5][1],t,p[4],p[7],view(du.x[2],1,:))

    Tel_func(u.x[3][2:end],u.x[4][2:end],p[2],p[2].DOS,p[3],p[1],p[5][2:end],t,p[6][2:end],p[4],view(du.x[3],2:p[4].length))
    Tph_func(u.x[3][2:end],u.x[4][2:end],p[2],p[2].DOS,p[3],p[5][2:end],view(du.x[4],2:p[4].length))
    du.x[1][2:end] .= 0.0

    nothing
end
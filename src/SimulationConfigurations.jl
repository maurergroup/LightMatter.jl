function run_simulation(key_list,initialtemps,tspan,sim,mp,las,dim,cons;save=2.0,tolerance=1e-4,min_step=0.1)
    u0 = generate_initialconditions(key_list,mp,initialtemps,dim)
    p = generate_parameters(sim,mp,cons,las,initialtemps,dim)
    return run_dynamics(p,u0,tspan,save,tolerance,min_step)
end

function run_dynamics(p,u0,tspan,save,tolerance,min_step)
    if length(u0.x) == 1
        return run_dynamics1(p,u0,tspan,save,tolerance,min_step)
    elseif length(u0.x) == 2
        return run_dynamics2(p,u0,tspan,save,tolerance,min_step)
    elseif length(u0.x) == 3
        return sol = run_dynamics3(p,u0,tspan,save,tolerance,min_step)
    elseif length(u0.x) == 4
        return sol = run_dynamics4(p,u0,tspan,save,tolerance,min_step)
    end
end

function run_dynamics4(p,u0,tspan,save,tolerance,min_step)
    println("Running small time span to precompile dynamics")
    solve(ODEProblem(FullAthEM_simulation,u0,(0.0,0.1),p),Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmin=min_step)
    println("Running main dynamics")
    prob=ODEProblem(FullAthEM_simulation,u0,tspan,p)
    sol = solve(prob,Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmin=min_step)
    return sol
end

function run_dynamics1(p,u0,tspan,save,tolerance,min_step)
    println("Running small time span to precompile dynamics")
    solve(ODEProblem(EHP_simulation,u0,(0.0,0.1),p),Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmin=min_step)
    println("Running main dynamics")
    prob=ODEProblem(EHP_simulation,u0,tspan,p)
    sol = solve(prob,Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmin=min_step)
    return sol
end
   
function run_dynamics3(p,u0,tspan,save,tolerance,min_step)
    println("Running small time span to precompile dynamics")
    solve(ODEProblem(eeAthEM_simulation,u0,(0.0,0.1),p),Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmin=min_step)
    println("Running main dynamics")
    prob=ODEProblem(eeAthEM_simulation,u0,tspan,p)
    sol = solve(prob,Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmin=min_step)
    return sol
end

function run_dynamics2(p,u0,tspan,save,tolerance,min_step)
    println("Running small time span to precompile dynamics")
    solve(ODEProblem(TTM_simulation,u0,(0.0,0.1),p),Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmin=min_step)
    println("Running main dynamics")
    prob=ODEProblem(TTM_simulation,u0,tspan,p)
    sol = solve(prob,Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmin=min_step)
    return sol
end

function TTM_simulation(du,u,p,t)
    println(t)
    @assert all(!isnan(x) for x in u)
    electrontemperature_conductivity(u.x[1],p[4],u.x[2],p[2],p[7])
    Threads.@threads for i in eachindex(u.x[1])
        @inbounds p[6][i] = find_chemicalpotential(p[5],u.x[1][i],p[2].DOS[i],p[3].kB,p[2].FE,p[2].n0)
    end
    Tel_func(u.x[1],u.x[2],p[2],p[2].DOS,p[3],p[1],p[6],t,p[7],p[4],du.x[1])
    Tph_func(u.x[1],u.x[2],p[2],p[2].DOS,p[3],p[6],du.x[2])
    nothing
end

function eeAthEM_simulation(du,u,p,t)
    println(t) 
    #n = 1, Tel = 2, fneq = 3 

    μ = find_chemicalpotential.(u.x[1],u.x[2],Ref(p[2].DOS[1]),p[3].kB,p[2].FE,p[2].n0)
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
        p[5][i] = find_chemicalpotential(u.x[2][i],u.x[1][i],p[2].DOS[i],p[3].kB,p[2].FE,p[2].n0)
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
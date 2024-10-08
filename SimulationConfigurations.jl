function run_simulation(key_list,initialtemps,tspan)
    sim,mp,las,dim,cons = setup()
    u0 = generate_initalconditions(key_list,mp,cons,initialtemps,dim)
    p = generate_parameters(sim,mp,cons,las,initialtemps,dim)
    return run_dynamics(p,u0,tspan)
end

function run_dynamics(p,u0,tspan)
    if length(u0.x) == 1
        return run_dynamics1(p,u0,tspan)
    elseif length(u0.x) == 2
        return run_dynamics2(p,u0,tspan)
    elseif length(u0.x) == 3
        return sol = run_dynamics3(p,u0,tspan)
    elseif length(u0.x) == 4
        return sol = run_dynamics4(p,u0,tspan)
    end
end

function run_dynamics4(p,u0,tspan)
    prob=ODEProblem(FullAthEM_simulation,u0,tspan,p)
    sol = solve(prob,Tsit5(),abstol=1e-4,reltol=1e-4,saveat=1.0,dtmin=0.1)
    return sol
end

function run_dynamics1(p,u0,tspan)
    prob=ODEProblem(EHP_simulation,u0,tspan,p)
    sol = solve(prob,Tsit5(),abstol=1e-4,reltol=1e-4,saveat=1.0,dtmin=0.1)
    return sol
end
   
function run_dynamics3(p,u0,tspan)
    prob=ODEProblem(AthEM_simulation,u0,tspan,p)
    sol = solve(prob,OrdinaryDiffEq.Tsit5(),abstol=1e-4,reltol=1e-4,saveat=1.0,dtmin=0.1)
    return sol
end

function run_dynamics2(p,u0,tspan)
    prob=ODEProblem(TTM_simulation,u0,tspan,p)
    sol = solve(prob,Tsit5(),abstol=1e-4,reltol=1e-4,saveat=2.0,dtmin=0.1)
    return sol
end

function TTM_simulation(du,u,p,t)
    println(t)
    electrontemperature_conductivity(u.x[1],p[4],u.x[2],p[2],p[7])
    Threads.@threads for i in eachindex(u.x[1])
        @inbounds p[6][i] = find_chemicalpotential(p[5],u.x[1][i],p[2].DOS,p[3].kB,p[2].FE,p[2].n0)
    end

    Tel(u.x[1],u.x[2],p[2],p[3],p[1],p[6],t,p[7],p[4],du.x[1])
    Tph(u.x[1],u.x[2],p[2],p[3],p[6],du.x[2])
    nothing
end

function AthEM_simulation(du,u,p,t)
    println(t) 
    #n = 1, Tel = 2, fneq = 3 

    μ = zeros(length(u.x[1]))
    Threads.@threads for i in eachindex(u.x[1])
        μ[i] = find_chemicalpotential(u.x[1][i],u.x[2][i],p[2].DOS,p[3].kB,p[2].FE,p[2].n0)
    end
    relax_dis = zeros(1,length(u.x[3]))
    relax(u.x[2],u.x[3],u.x[1],μ,p[2],p[3],relax_dis)
    noe(relax_dis,μ,p[2],du.x[1])
    Tel(u.x[2],p[2],p[3],μ,relax_dis,du.x[1],zeros(length(u.x[2])),du.x[2])
    fneq(u.x[3],u.x[2],p[2],p[3],p[1],μ,t,p[4],relax_dis,du.x[3])
    nothing
end

function FullAthEM_simulation(du,u,p,t)
    println(t)

    electrontemperature_conductivity(u.x[1],p[4],u.x[3],p[2],p[6])
    Threads.@threads for i in eachindex(p[5])
        p[5][i] = find_chemicalpotential(u.x[2][i],u.x[1][i],p[2].DOS,p[3].kB,p[2].FE,p[2].n0)
    end
    relax(u.x[1],u.x[4],u.x[2],p[5],p[2],p[3],p[7])
    noe(p[7],p[5],p[2],du.x[2])
    Tel(u.x[1],u.x[3],p[2],p[3],p[5],p[7],du.x[2],p[6],du.x[1])
    Tph(u.x[1],u.x[3],p[2],p[3],p[5],u.x[4],du.x[3])
    fneq(u.x[4],u.x[1],p[2],p[3],p[1],p[5],t,p[4],p[7],du.x[4])
    nothing
end

function EHP_simulation(du,u,p,t)
    println(t)
    vars = (Tel=p[7],fneq=u.x[neq_idx][:],n=p[8],μ=p[9],cons=p[5],mp=p[4],t=t,las=p[3])
    du.x[neq_idx][:].=eval_expr(p[1]["fneq"],vars)
end


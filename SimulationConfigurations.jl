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
    sol = solve(prob,Tsit5(),abstol=1e-4,reltol=1e-4,saveat=1.0,dtmin=0.1)
    return sol
end

function run_dynamics2(p,u0,tspan)
    prob=ODEProblem(TTM_simulation,u0,tspan,p)
    sol = solve(prob,Tsit5(),abstol=1e-4,reltol=1e-4,saveat=1.0,dtmin=0.1)
    return sol
end

function TTM_simulation(du,u,p,t)
    println(t)
    
    conducted_Tel=electrontemperature_conductivity(u.x[1],p[4],u.x[2],p[2])
    μ = find_chemicalpotential.(p[5],u.x[1],Ref(p[2].DOS),p[3].kB,p[2].FE)

    Tel(u.x[1],u.x[2],p[2],p[3],p[1],μ,t,conducted_Tel,du.x[1])
    Tph(u.x[1],u.x[2],p[2],p[3],μ,du.x[2])
end

function AthEM_simulation(du,u,p,t)
    println(t)
    Tel_idx = findfirst(==("Tel"),p[2])
    neq_idx = findfirst(==("fneq"),p[2])
    n_idx = findfirst(==("noe"),p[2])

    μ = find_chemicalpotential(u.x[n_idx][1],u.x[Tel_idx][1],p[4].DOS,p[5].kB,p[4].FE)

    relax_vars = (Tel=u.x[Tel_idx][1],fneq=u.x[neq_idx][:],n=u.x[n_idx][1],μ=μ,cons=p[5],mp=p[4],t=t,las=p[3])
    relax_dis = eval_expr(athem_electronelectronscattering(),relax_vars)

    n_vars = merge(relax_vars,(Symbol("relax_dis")=>relax_dis,))
    du.x[n_idx][:] .= eval_expr(p[1]["noe"],n_vars)
    vars = merge(n_vars,(Symbol("Δn")=>du.x[n_idx][1][1],))
    du.x[Tel_idx][:].=eval_expr(p[1]["Tel"],vars)
    du.x[neq_idx][:] .= eval_expr(p[1]["fneq"],vars)
end

function FullAthEM_simulation(du,u,p,t)
    println(t)
    Tel_idx = findfirst(==("Tel"),p[2])
    Tph_idx = findfirst(==("Tph"),p[2])
    neq_idx = findfirst(==("fneq"),p[2])
    n_idx = findfirst(==("noe"),p[2])

    conducted_Tel=electrontemperature_conductivity(u.x[Tel_idx],p[6],u.x[Tph_idx],p[4])

    Threads.@threads for i in eachindex(u.x[Tel_idx])
        μ = find_chemicalpotential(u.x[n_idx][i],u.x[Tel_idx][i],p[4].DOS,p[5].kB,p[4].FE)

        relax_vars = (Tel=u.x[Tel_idx][i],fneq=u.x[neq_idx][i,:],n=u.x[n_idx][i],μ=μ,cons=p[5],mp=p[4]
        ,t=t,las=p[3],Tph=u.x[Tph_idx][i],cond=conducted_Tel[i])
        relax_dis = eval_expr(athem_electronelectronscattering(),relax_vars)

        n_vars = merge(relax_vars,(Symbol("relax_dis")=>relax_dis,))
        du.x[n_idx][i] .= eval_expr(p[1]["noe"],n_vars)

        vars = merge(n_vars,(Symbol("Δn")=>du.x[n_idx][1][1],))
        du.x[Tel_idx][i].=eval_expr(p[1]["Tel"],vars)
        du.x[neq_idx][i,:].=eval_expr(p[1]["fneq"][i],vars)
        du.x[Tph_idx][i].=eval_expr(p[1]["Tph"],vars)
    end
end

function EHP_simulation(du,u,p,t)
    println(t)
    vars = (Tel=p[7],fneq=u.x[neq_idx][:],n=p[8],μ=p[9],cons=p[5],mp=p[4],t=t,las=p[3])
    du.x[neq_idx][:].=eval_expr(p[1]["fneq"],vars)
end


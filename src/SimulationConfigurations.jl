function run_simulation(sys,initialtemps,tspan,sim,mp,las,dim;save=2.0,tolerance=1e-4,max_step=0.1,min_step=0.01)
    u0 = generate_initialconditions(sim,mp,initialtemps,dim)
    p = generate_parameters(sim,las,mp,initialtemps,dim)
    simulation_expr = simulation_construction(sys,sim)
    simulation_problem = mk_function((:du,:u,:p,:t),(),simulation_expr)
    println("Running small time span to precompile dynamics")
    solve(ODEProblem(simulation_problem,u0,(0.0,0.1),p),Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmax=max_step,dtmin=min_step)
    println("Running main dynamics")
    prob=ODEProblem(simulation_problem,u0,tspan,p)
    sol = solve(prob,Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmax=max_step,dtmin=min_step)
    return sol
end

function run_dynamics(p,u0,tspan,save,tolerance,max_step,min_step)
    println("Running small time span to precompile dynamics")
    solve(ODEProblem(simulation_problem,u0,(0.0,0.1),p),Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmax=max_step,dtmin=min_step)
    println("Running main dynamics")
    prob=ODEProblem(simulation_problem,u0,tspan,p)
    sol = solve(prob,Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmax=max_step,dtmin=min_step)
    return sol
end
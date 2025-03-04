"""
    run_simulation(sys::Dict{String,Union{Expr,Vector{Expr}}},initialtemps::Dict{String, Real},
    tspan::Tuple{Real,Real},sim::SimulationSettings,mp::MaterialParameters,las::Laser,dim::Dimension
    ;save=2.0,tolerance=1e-4,max_step=0.1,min_step=0.01)

    Creates the ODE simulation that will be performed as well as performs the ODE itself returning a DiffEq.jl
    solution object with the final results inside. It first runs a small simulation between 0 and 0.1 to precompile
    all functions before running the full simulation. Currently, it uses an adaptable time-stepping Runge-Kutta method to 
    perform the dynamics but this may be flexible in the future. Tolerance defines both the absolute and relative tolerance. 
"""
function run_simulation(sys::Dict{String,Union{Expr,Vector{Expr}}},initialtemps::Dict{String, <:Real},
    tspan::Tuple{Real,Real},sim::SimulationSettings,mp::MaterialParameters,las::Laser,dim::Dimension
    ;save=2.0,tolerance=1e-4,max_step=0.1,min_step=0.01)

    u0 = generate_initialconditions(sim,mp,initialtemps,dim)
    p = generate_parameters(sim,las,mp,initialtemps,dim)
    simulation_expr = simulation_construction(sys,sim)
    simulation_problem! = mk_function((:du,:u,:p,:t),(),simulation_expr)
    println("Running small time span to precompile dynamics")
    solve(ODEProblem(simulation_problem!,u0,(0.0,0.1),p),Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmax=max_step,dtmin=min_step)
    println("Running main dynamics")
    prob=ODEProblem(simulation_problem!,u0,tspan,p)
    sol = solve(prob,Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmax=max_step,dtmin=min_step)
    return sol
end
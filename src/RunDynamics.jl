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
    tspan::Tuple{Real,Real},sim::Simulation
    ;save=2.0,tolerance=1e-4,max_step=0.1,min_step=0.01,callbacks=CallbackSet())

    u0 = generate_initialconditions(sim,initialtemps)
    p = generate_parameters(sim,initialtemps)
    simulation_expr = simulation_construction(sys,sim)
    simulation_problem! = mk_function((:du,:u,:p,:t),(),simulation_expr)
    println("Precompiling")
    simulation_problem!(similar(u0),u0,p,0.0)
    println("Running main dynamics")
    prob=ODEProblem(simulation_problem!,u0,tspan,p)
    sol = solve(prob,Tsit5(),abstol=tolerance,reltol=tolerance,saveat=save,dtmax=max_step,dtmin=min_step,callback = callbacks)
    return sol
end
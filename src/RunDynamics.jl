###
# Could add some more key-word arguemnts such as the algorithm used
###
"""
    run_simulation(sys::Dict{String,Union{Expr,Vector{Expr}}}, initialtemps::Dict{String, <:Real},
    tspan::Tuple{Real,Real}, sim::Simulation; 
    save, tolerance, max_step, min_step, callbacks)
    
    Generates the problem the dynamics will solve and then solves the coupled system of ODE's.
    Currently always uses Tsit5 for the integration routine but in the future that may be user-defined

    # Arguments
    - 'sys': Dictionary of ODE equations to be propagated
    - 'initialtemps': Dictionary of initial temperatures of the bath
    - 'tspan': Tuple of values for the dynamics to run between (the laser is centred on 0.0)
    - 'sim': Simulation settings and parameters
    
    # Optional Arguments
    - 'save': The time points the solution is saved at : Defaults to 2.0 fs
    - 'tolerance': The absolute and relative tolerance value of the dynamics : Defaults to 1e-4
    - 'max_step': Maximum step size for adaptive integrator : Defaults to 0.1 fs
    - 'min_step': Minimum step size for adaptive integrator : Defaults to 0.01 fs
    = 'callbacks': Any user-defined callbacks : Defaults to CallbackSet()

    # Returns
    - The solution of the dynamics calculation
"""
function run_simulation(sys::Dict{String,Union{Expr,Vector{Expr}}}, initialtemps::Dict{String, <:Real},
    tspan::Tuple{Real,Real}, sim::Simulation; 
    save=2.0, tolerance=1e-4, max_step=0.1, min_step=0.01, callbacks=CallbackSet())

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
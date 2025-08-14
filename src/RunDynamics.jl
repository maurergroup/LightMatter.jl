"""
    run_simulation(sys::Dict{String,Union{Expr,Vector{Expr}}}, initialtemps::Dict{String, Float64},
    tspan::Tuple{Float64,Float64}, sim::Simulation; 
    save, tolerance, max_step, min_step, callbacks)
    
    Generates the problem the dynamics will solve and then solves the coupled system of ODE's.
    Currently always uses Tsit5 for the integration routine but in the future that may be user-defined

    # Arguments
    - 'sys': Dictionary of ODE equations to be propagated
    - 'initialtemps': Dictionary of initial temperatures of the bath
    - 'tspan': Tuple of values for the dynamics to run between (the laser is centred on 0.0)
    - 'sim': Simulation settings and parameters
    
    # KWARGS
    - Any key-word arguemnts from DiffEq that work with ODEProblems can be including in a namedtuple here

    # Returns
    - The solution of the dynamics calculation
"""
function run_simulation(sim::Simulation, initialtemps::Dict{String, Float64},
    tspan::Tuple{Float64,Float64}; alg=Tsit5(), print_time=false, kwargs...)
    
    sys = function_builder(sim)
    u0 = generate_initialconditions(sim,initialtemps)
    p = generate_parameters(sim,initialtemps)
    simulation_expr = simulation_construction(sys,sim, print_time)
    simulation_problem! = mk_function((:du,:u,:p,:t),(),simulation_expr)
    println("Precompiling") 
    simulation_problem!(similar(u0),u0,p,0.0)
    prob=ODEProblem(simulation_problem!,u0,tspan,p)
    println("Running Script")
    sol = solve(prob, alg; kwargs...)
    return sol
end

# Temporary function overloading to work with stiff integrators such as Trapezoid()

using ArrayInterface

function ArrayInterface.zeromatrix(A::NamedArrayPartition)
   B = ArrayPartition(A)
    x = reduce(vcat,vec.(B.x))
    x .* x' .* false
end
using ModelingToolkit,DifferentialEquations
using ModelingToolkit: t_nounits as t, D_nounits as D

function eq_factory(length;name)
    @variables y(t)[1:length]

    @named dy = y_eq_template(length)

    connections = [dy.rhs ~ rhs_function_wrapper(dy.y,length),
                   dy.y ~ y]
    
    compose(ODESystem(connections,t;name),dy)
end

function y_eq_template(length;name)
    @variables y(t)[1:length] rhs(t)[1:length]

    eqs = D(y) ~ rhs

    ODESystem(eqs,t;name)
end

function rhs_function_wrapper(y,length)
    @parameters a b

    return rhs_function(y,a,b,t)
end

function rhs_function(y::AbstractVector,a::Real,b::Real,t::Real)
    return ((y*a).+b)./(t)
end
@register_array_symbolic rhs_function(y::AbstractVector,a::Num,b::Num,t::Num) begin
    size = (length(y),)
    eltype = eltype(y)
end

@named test_sys = eq_factory(10)
test_simp = structural_simplify(test_sys)
#= u0 = [test_simp.y=>zeros(10)]
p = [test_simp.a=>2,
     test_simp.b=>5]
tspan=(0.0,10.0)
prob=ODEProblem(test_simp,u0,tspan,p) =#
#sol = solve(prob,Rosenbrock23();reltol=1e-3,abstol=1e-3)


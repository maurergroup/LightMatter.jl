using ModelingToolkit,DifferentialEquations,Plots,Symbolics,Dierckx,DelimitedFiles,Integrals
using BenchmarkTools,ForwardDiff,StaticArrays,IfElse,Cubature, SimpleDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

function eq_factory(source::Num,egl::Int;name)
    @parameters (xgrid)[1:egl] T a
    @named Deq = eq_template(egl)
    
    connections = [Deq.rhs ~ symbolic_func(xgrid,T,Deq.y,a).*source]

    compose(ODESystem(connections,t;name),Deq)
end

function eq_template(egl::Int;name)
    @variables (y(t))[1:egl] (rhs(t))[1:egl]

    eqs = D(y) ~ rhs

    ODESystem(eqs,t;name)
end

function symbolic_func(xgrid,T,y,a)

    tot = y.+exp.(xgrid./T)
    totspl=get_interpolate(xgrid,tot)

    Δplus = delta_plus(xgrid,totspl,a)
    Δminus = delta_minus(xgrid,totspl,a)

    pc = delta_integral(Δplus,xgrid) / delta_integral(Δminus,xgrid)
    Δshape = (Δplus.*pc).-Δminus
    inten = x_integral(Δshape,xgrid)

    return Δshape./inten
end
@register_array_symbolic  symbolic_func(xgrid::AbstractVector,T::Num,y::AbstractVector,a::Num) begin
    size = (length(xgrid),)
    eltype = eltype(y)
end

function delta_plus(xgrid,totspl,a)
    return totspl(xgrid.-a).*(1 .-totspl(xgrid))
end

function delta_minus(xgrid,totspl,a)
    return totspl(xgrid).*(1 .-totspl(xgrid.+a))
end

get_interpolate(xvals::AbstractVector,yvals::AbstractVector) = Spline1D(xvals,yvals,bc="nearest")

function x_integral(Dis,xgrid)
    integrand = Dis.*xgrid
    prob = SampledIntegralProblem(integrand,xgrid)
    return solve(prob,SimpsonsRule()).u
end

function delta_integral(Dis::AbstractVector,xgrid::AbstractVector)::Real
    integrand = Dis
    prob = SampledIntegralProblem(integrand,xgrid)
    return solve(prob,SimpsonsRule()).u
end

function eulermethod(time,xgrid,T,y,a,dt,source)
    input = substitute(source,Dict([h=>62.0,σ=>50.0,t=>time]))
    Δ = symbolic_func(xgrid,T,y,a)
    return Δ*dt*input.val
end

@parameters h σ 
ex = h/(σ*sqrt(2*pi))*exp(-t^2/σ^2)
grid = range(-6,6,step=0.1)
egl=length(grid)

y = zeros(egl)
@btime for time in range(-250,250,step=0.1)
    y .+= eulermethod(time,grid,0.3,y,3.0,0.1,ex)
end

@named test_eq = eq_factory(ex,egl)
simp = structural_simplify(test_eq)

u0=[simp.Deq.y=>zeros(egl)]
p=[simp.xgrid=>grid,
simp.h=>62.0,
simp.σ=>50.0,
simp.T=>0.3,
simp.a=>3.0]
prob=ODEProblem(simp,u0,(-250.0,250.0),p)
@btime sol=solve(prob,SimpleATsit5(),dt=0.1;abstol=1e-3,reltol=1e-3)
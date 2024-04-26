module ETTM
include("StockEquations.jl")
using Interpolations,Integrals,Roots,.Equations

function plasma(mp,gc) ::Float64
    mp.ω= sqrt(mp.n*mp.ne*gc.eleccharge^2/(mp.effmass*mp.elecmass*gc.ϵ0))
end

function tau_0(mp)::Float64
    sv.τ= 128/(sqrt(3)*pi^2*mp.ω)
end

function tau_ee(E::Float64,sv,gc,mp)::Float64
    return sv.τ*sv.μ^2 /((E-sv.μ)^2 +(pi*gc.kB*sv.Tel)^2)
end

function tau_ep()::Float64
    τep=τf*hv/(kB*θ)
    return τep
end

function tau_th(mp,sv)::Float64
    return (mp.γ*(sv.Tel+sv.Tph)/(2*mp.g))
end

function probint(sv,E_h::Float64,E_e::Float64,gc)::Float64
    return sv.DOS(E_h)*sv.DOS(E_e)*Equations.FermiDirac(E_h,sv,gc)*(1-Equations.FermiDirac(E_e,sv,gc))#*gaussian(E_e-E_h)
end

function particlenumber(vec::Vector{Float64},sv,ERange::Vector{Float64})::Float64
    return sum(vec.*sv.DOS.(ERange).*(ERange[2]-ERange[1]))
end

function particleconstant(sv,l,gc)::Float64
    ERange=range(-4+sv.μ,4+sv.μ,step=0.01)
    elec=probint.(Ref(sv),ERange.-l.hv,ERange,Ref(gc))
    hol=probint.(Ref(sv),ERange,ERange.+l.hv,Ref(gc))
    f(x)=(particlenumber(elec*x,sv,ERange))-particlenumber(hol,sv,ERange)
    x0=1
    prob=ZeroProblem(f,x0)
    return solve(prob,Order8();atol=1e-3,rtol=1e-3)
end

function FDChange(E::Float64,sv,particle::Float64,gc)::Float64
    elec=probint(sv,E-hv,E,gc)
    hol=probint(sb,E,E+hv,gc)
    return (particle*elec).-hol
end

function stepint(sv,gc,particle::Float64)::Float64
    int(u,p)=FDChange(u,sv,particle,gc)*sv.DOS(u)*u
    prob=IntegralProblem(int,sv.μ-(3*hv),sv.μ+(3*hv))
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-7,abstol=1e-7)
    return sol.u
end

function stepsize(sv,gc,tprime::Float64,particle::Float64,l)::Float64
    return Equations.LaserStrength(tprime,l,mp)/stepint(sv,gc,particle)
end

function ueeintegral!(u::Float64,sv,t::Float64,tprime::Float64,gc,mp,particle::Float64)::Float64
    #1/τee*FDChange*exp(-((t-tprime)/τee)-((t-tprime)/τep))*DOS(E)*E
    τee=tau_ee(E::Float64,sv,gc,mp)
    return τee*FDChange(u,sv,particle,gc)*exp((t-tprime)/τee)*sv.DOS(u)*u
end

function instantuee(sv,mp,l,gc,tprime::Float64,t::Float64,particle::Float64,τep::Float64)::Float64
    int(u,p)=ueeintegral!(u,sv,t,tprime,gc,mp,particle)
    prob=IntegralProblem(int,sv.μ-(3*l.hv),sv.μ+(3*l.hv))
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-5,abstol=1e-5)
    return abs(stepsize(sv,gc,tprime,particle,l))*sol.u
end

function totalUee(sv,gc,l,mp,particle::Float64,τep::Float64,t::Float64)::Float64
    int(u,p)=instantuee(sv,mp,l,gc,u,t,particle,τep)
    prob=IntegralProblem(int,0.0,t)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-3,abstol=1e-3)
    return sol.u
end

function uepintegral!(u::Float64,sv,t::Float64,tprime::Float64,gc,mp,particle::Float64,τep::Float64)::Float64
    #1/τee*FDChange*exp(-((t-tprime)/τee)-((t-tprime)/τep))*DOS(E)*E
    τee=tau_ee(E::Float64,sv,gc,mp)
    return 1/τep*FDChange(u,sv,particle,gc)*exp(((t-tprime)/τee)-((t-tprime)/τep))*sv.DOS(u)*u
end

function instantuep(sv,mp,l,gc,tprime::Float64,t::Float64,particle::Float64,τep)::Float64
    int(u,p)=uepintegral!(u,sv,t,tprime,gc,mp,particle,τep)
    prob=IntegralProblem(int,sv.μ-(3*l.hv),sv.μ+(3*l.hv))
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-5,abstol=1e-5)
    return abs(stepsize(sv,gc,tprime,particle,l))*sol.u
end

function totalUep(sv,gc,l,mp,particle::Float64,τep::Float64,t::Float64)::Float64
    int(u,p)=instantuep(sv,mp,l,gc,u,t,particle,τep)
    prob=IntegralProblem(int,0.0,t)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-3,abstol=1e-3)
    return sol.u
end

export particleconstant

end
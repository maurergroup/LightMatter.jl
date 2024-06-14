using ModelingToolkit,DifferentialEquations,Plots,Symbolics,Dierckx,DelimitedFiles,Integrals
using NonlinearSolve,ForwardDiff,QuadGK,StaticArrays,BenchmarkTools,Unitful
using ModelingToolkit: t_nounits as t, D_nounits as D 

include("../SymbolicsInterpolation.jl")

abstract type Laser end
@kwdef struct Gaussian <: Laser
    FWHM::Real
    Offset::Real
    Power::Real
    hv::Real
    Reflectivity::Real
    Transport::String
end

mutable struct params
    g::Float64
    Cph::Float64
    ϵ::Float64
    kB::Float64
    n::Float64
    ϕ::Float64
    Offset::Float64
    FWHM::Float64
end

get_interpolate(xvals::Vector{Float64},yvals::Vector{Float64}) = Spline1D(xvals,yvals,bc="nearest")


function generate_DOS(File::String,FE,n)
    TotalDOS::Matrix{Float64}=readdlm(File,skipstart=3)
    return get_interpolate(TotalDOS[:,1].+FE,TotalDOS[:,2].*n)
end

function (::Gaussian)()
    @parameters FWHM Offset
    temp = sqrt(4*log(2)/pi)/FWHM*exp(-4*log(2)*((t-(2*FWHM)-Offset)/FWHM)^2)
    return temp
end

function dTel(;name)
    @variables S(t) Tel(t) HC(t) EP(t)

    eqs = D(Tel) ~ (-EP + S)/HC

    ODESystem(eqs,t;name)

end

function dTph(;name)
    @parameters g Cph
    @variables Tel(t) Tph(t)
    
    eqs = D(Tph) ~ g*(Tel-Tph)/Cph

    ODESystem(eqs,t;name)

end

function Laser(lp::Gaussian)
    @parameters ϵ ϕ
    return  lp()*ϕ/ϵ
end

function dTel_factory(DOS::Spline1D,lp::Gaussian;name)
    @parameters kB μ g
    @variables Tph(t)
    @named dTeldt = dTel()
    laser=Laser(lp)
    connections=[dTeldt.S ~ laser,
                 dTeldt.HC ~ HeatCapacity(kB,dTeldt.Tel,μ,DOS),
                 dTeldt.EP ~ g*(dTeldt.Tel-Tph)]
    compose(ODESystem(connections,t;name),dTeldt)
end

function dFDdT(kB::Float64,Tel::Float64,μ::Float64,E::Float64)::Real
    Numer=(E-μ)*exp((E-μ)/(kB*Tel))
    Denom=kB*Tel^2*(exp((E-μ)/(kB*Tel))+1)^2
    return Numer/Denom
end

function HeatCapacity_int(u::Float64,p::Tuple{Float64,Float64,Float64,Spline1D})
    return abs(dFDdT(p[1],p[2],p[3],u)*p[4](u)*u)
end

function HeatCapacity(kB::Float64,Tel::Float64,μ::Float64,DOS::Spline1D)
    p=(kB,Tel,μ,DOS)
    int(u,p) = HeatCapacity_int(u,p)
    return solve(IntegralProblem(int,(μ-10,μ+10),p),HCubatureJL(initdiv=10);reltol=1e-5,abstol=1e-5).u
end
@register_symbolic HeatCapacity(kB::Num,Tel::Num,μ::Num,DOS::Spline1D)

function get_FermiEnergy(File)
    TotalDOS::Matrix{Float64}=readdlm(File,skipstart=3)
    nonzero=findfirst(TotalDOS[:,2] .== filter(x -> x!=0,TotalDOS[:,2])[1])
    return abs(TotalDOS[nonzero,1])
end

function get_noparticles_int(y::Vector{Float64},u::Vector{Float64},p::Tuple{Float64,Float64,Float64,Spline1D})
    n=Threads.nthreads()
    Threads.@threads for i in 1:n
        y[i:n:end] .= FermiDirac.(@view(u[i:n:end]),p[1],p[2],p[3]).*p[4](@view(u[i:n:end]))
    end
end

function get_noparticles_int(y::Vector{ForwardDiff.Dual},u::Vector{Float64},p::Tuple{ForwardDiff.Dual,Float64,Float64,Spline1D})
    n=Threads.nthreads()
    Threads.@threads for i in 1:n
        y[i:n:end] .= FermiDirac.(@view(u[i:n:end]),p[1],p[2],p[3]).*p[4](@view(u[i:n:end]))
    end
end

function get_noparticles(μ::Float64,Tel::Float64,DOS::Spline1D,kB::Float64)
    p=(μ,Tel,kB,DOS)
    int = BatchIntegralFunction(get_noparticles_int,zeros(0))
    return solve(IntegralProblem(int,(μ-10,μ+10),p),QuadGKJL();abstol=1e-3).u
end

function get_noparticles(μ::ForwardDiff.Dual,Tel::Float64,DOS::Spline1D,kB::Float64)
    p=(μ,Tel,kB,DOS)
    int = BatchIntegralFunction(get_noparticles_int,Array{ForwardDiff.Dual}(undef,0))
    return solve(IntegralProblem(int,(ForwardDiff.value(μ)-10,ForwardDiff.value(μ)+10),p),QuadGKJL();reltol=1e-3,abstol=1e-3).u
end

FermiDirac(E::Float64,μ::Float64,Tel::Float64,kB::Float64) = 1/(exp((E-μ)/(kB*Tel))+1)

FermiDirac(E::Float64,μ::ForwardDiff.Dual,Tel::Float64,kB::Float64) = 1/(exp((E-μ)/(kB*Tel))+1)


function find_chemicalpotential(no_part::Float64,Tel::Float64,μ::Float64,DOS::Spline1D,kB::Float64)
    f(u,p) = no_part - get_noparticles(u,Tel,DOS,kB)
    sol = solve(NonlinearProblem(f,μ),SimpleKlement();abstol=1e-3,reltol=1e-3)
    return sol.u
end
@register_symbolic find_chemicalpotential(no_part::Num,Tel::Num,μ::Num,DOS::Spline1D,kB::Num)

function update_chempot!(integ,u,p,ctx::Tuple{Spline1D,Float64})
    integ.p[p.μ] = find_chemicalpotential(ctx[2],integ.u[u.Tel],integ.p[p.μ],ctx[1],integ.p[p.kB])
end

function unitconversion(stct::params)
    stct.g = ustrip(uconvert(u"eV/nm^3/fs/K",stct.g*u"W/m^3/K"))
    stct.Cph = ustrip(uconvert(u"eV/nm^3/K",stct.Cph*u"J/m^3/K"))
    stct.ϵ = ustrip(uconvert(u"nm",stct.ϵ*u"m"))
    stct.kB = ustrip(uconvert(u"eV/K",stct.kB*u"J/K"))
    stct.n = ustrip(uconvert(u"nm^-3",stct.n*u"m^-3"))
    stct.ϕ = ustrip(uconvert(u"eV/nm^2",stct.ϕ*u"J/m^2"))
    stct.Offset = ustrip(uconvert(u"fs",stct.Offset*u"s"))
    stct.FWHM = ustrip(uconvert(u"fs",stct.FWHM*u"s"))
    return stct
end

function main()
    
    param = params(2.3e16,2.5e6,1.27e-8,1.380649e-23,5.9e28,10.0,200e-15,50e-15)
    param = unitconversion(param)
    param.Offset = param.Offset-(2*param.FWHM)
    DOS=generate_DOS("DOS/Au_DOS.dat",0.0,param.n)
    lp = Gaussian(FWHM=param.FWHM,Offset=param.Offset,Power=param.ϕ,hv=3.1,Reflectivity=0.0,
    Transport="Optical")
    param_storage=[]

    @named Tph_eq = dTph() 
    @named Tel_eq=dTel_factory(DOS,lp)
    no_part = get_noparticles(0.0,1e-16,DOS,param.kB)
    chempot = (t>=0.0) => (update_chempot!,[Tel_eq.dTeldt.Tel=>:Tel],
    [Tel_eq.μ=>:μ,Tel_eq.kB=>:kB],[Tel_eq.μ],(DOS,no_part))
    connections=[Tel_eq.Tph ~ Tph_eq.Tph,
                Tph_eq.Tel ~ Tel_eq.dTeldt.Tel]
    connected = compose(ODESystem(connections,t,name=:connected
                ,defaults=Pair{Num,Any}[Tel_eq.g => Tph_eq.g],discrete_events=chempot),Tel_eq,Tph_eq)
    connected_simp=structural_simplify(connected)

    u0=[Tel_eq.dTeldt.Tel => 300.0,
        Tph_eq.Tph => 300.0]
    p=[Tph_eq.g => param.g,
      Tph_eq.Cph => param.Cph,
      Tel_eq.ϵ => param.ϵ,
      Tel_eq.kB => param.kB,
      Tel_eq.ϕ => param.ϕ,
      Tel_eq.Offset => param.Offset,
      Tel_eq.FWHM => param.FWHM,
      Tel_eq.μ => 0.0]

    prob=ODEProblem(connected_simp,u0,(0.0,400.0),p)
    sol=solve(prob,Tsit5(),abstol=1e-2,reltol=1e-2)
    plot(sol,idxs=[Tel_eq.dTeldt.EP,Tel_eq.dTeldt.S])
end

main()

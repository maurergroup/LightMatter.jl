using ModelingToolkit,DifferentialEquations,Plots,Symbolics,Dierckx,.SymbolicsInterpolation,DelimitedFiles,Integrals
using Unitful,BenchmarkTools,ForwardDiff
using ModelingToolkit: t_nounits as t, D_nounits as D 

include("SimulationSetup.jl")

function elec_template(;name)
    @variables Tel(t) Source(t) ElecPhon(t) HeatCapacity(t) Spatial(t)
    ParentScope(Tel)
    eqs = D(Tel) ~ (Source .+ Spatial .+ ElecPhon)./HeatCapacity

    ODESystem(eqs,t;name)
end

function t_electron_factory(mp::MaterialParameters,sim::SimulationSettings,laser::Num) 
    @named dTel = elec_template()
    connections =[dTel.HeatCapacity ~ t_electron_heatcapacity(mp,sim),
                  dTel.ElecPhon ~ t_electron_phononcoupling(mp,sim),
                  dTel.Spatial ~ 0.0,
                  dTel.Source ~ t_electron_sourceterm(laser,mp,sim)]
end

function t_electron_heatcapacity(mp::MaterialParameters,sim::SimulationSettings)
    if sim.ParameterApprox.ElectronHeatCapacity == true
        @parameters kB μ
        @variables Tel
        ParentScope(Tel)
        ParentScope(kB)
        ParentScope(μ)
        return nonlinear_electronheatcapacity(kB,Tel,μ,mp.DOS)
    else
        @parameters γ
        @variables Tel 
        ParentScope(γ)
        ParentScope(Tel)
        return γ*Tel
    end
end

function elec_heatcapacity_int(u::Float64,p::Tuple{Float64,Float64,Float64,Spline1D})
    return abs(dFDdT(p[1],p[2],p[3],u)*p[4](u)*u)
end

function nonlinear_elec_heatcapacity(kB::Float64,Tel::Float64,μ::Float64,DOS::Spline1D)
    p=(kB,Tel,μ,DOS)
    int(u,p) = HeatCapacity_int(u,p)
    return solve(IntegralProblem(int,(μ-(6*Tel/10000),μ+(6*Tel/10000)),p),HCubatureJL(initdiv=10);reltol=1e-5,abstol=1e-5).u
end
@register_symbolic nonlinear_elec_heatcapacity(kB::Num,Tel::Num,μ::Num,DOS::Spline1D)

function t_electron_phononcoupling(mp::MaterialParameters,sim::SimulationSettings)
    if sim.Interactions.ElectronPhonon == true
        if sim.ParameterApprox.ElectronPhononCoupling==true
            @parameters hbar kB λ μ
            @variables Tel Tph
            ParentScope(hbar)
            ParentScope(kB)
            ParentScope(λ)
            ParentScope(μ)
            ParentScope(Tel)
            ParentScope(Tph)
            return nonlinear_electronphononcoupling(hbar,kB,λ,mp.DOS,Tel,μ,Tph)
        else
            @parameters g
            @variables Tel Tph
            ParentScope(g)
            ParentScope(Tel)
            ParentScope(Tph)
            return -g*(Tel-Tph)
        end
    else
        return 0.0
    end
end

function electronphononcoupling_int(u::Float64,p::Tuple{Float64,Float64,Float64,Spline1D})
    return (p[4](u))^2*-1*dFDdE(p[1],p[2],p[3],u)
end

function electronphononcoupling_int(y::Vector{Float64},u::Vector{Float64},p::Tuple{Float64,Float64,Float64,Spline1D})
    n=Threads.nthreads()
    Threads.@threads for i in 1:n
        @inbounds y[i:n:end] .= p[4].(@view(u[i:n:end])).^2 .*-1 .*dFDdE.(p[1],p[2],p[3],@view(u[i:n:end]))
    end
end

function nonlinear_electronphononcoupling(hbar::Float64,kB::Float64,λ::Float64,DOS::Spline1D,Tel::Float64,μ::Float64,Tph::Float64)
    prefac=pi*kB*λ/DOS(μ)/hbar
    p=(kB,Tel,μ,DOS)
    int = BatchIntegralFunction(electronphononcoupling_int,zeros(0))
    g=prefac.*solve(IntegralProblem(int,(μ+(6*Tel/10000),μ+(6*Tel/10000)),p),QuadGKJL();reltol=1e-3,abstol=1e-3).u
    return g*(Tel-Tph)
end
@register_symbolic nonlinear_electronphononcoupling(hbar::Num,kB::Num,λ::Num,DOS::Spline1D,Tel::Num,μ::Num,Tph::Num)

function t_electron_sourceterm(laser::Num,mp::MaterialParameters,sim::SimulationSettings)
    if sim.Interactions.Nonequilibrium_electrons == true
        return electron_noneqelectron_coupling()
    else
        return laser
    end
end

function electron_noneqelectron_coupling()
end


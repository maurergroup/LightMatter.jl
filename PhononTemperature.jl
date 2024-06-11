using ModelingToolkit,DifferentialEquations,Plots,Symbolics,Dierckx,DelimitedFiles,Integrals
using Unitful,BenchmarkTools,ForwardDiff,StaticArrays
using ModelingToolkit: t_nounits as t, D_nounits as D 

include("SymbolicsInterpolation.jl")
include("SimulationVariables.jl")
include("SimulationSetup.jl")
include("Lasers.jl")

function t_phonon_factory(mp::MaterialParameters,sim::SimulationSettings;name)
    @variables Tph(t) Tel(t)
    @named dTph = t_phonon_template()
    connections=[dTph.Source ~ t_phonon_sourceterm(sim),
                 dTph.ElecPhon ~ -1*t_electron_phononcoupling(mp,sim),
                 dTph.HeatCapacity ~ t_phonon_heatcapacity(mp,sim),
                 dTph.Tph ~ Tph]
    compose(ODESystem(connections,t;name),dTph)
end

function t_phonon_template(;name)
    @variables Tph(t) Source(t) ElecPhon(t) HeatCapacity(t) Spatial(t)

    eqs = D(Tph) ~ (Source .+ ElecPhon)./HeatCapacity

    ODESystem(eqs,t;name)
end

function t_phonon_sourceterm(sim::SimulationSettings)
    if sim.Systems.NonEqElectrons == true
        if sim.Interactions.ElectronPhonon == true
            @variables uep(t)
            return uep
        else
            return 0.0
        end
    else
        return 0.0
    end
end

function t_phonon_heatcapacity(mp::MaterialParameters,sim::SimulationSettings)
    if sim.ParameterApprox.PhononHeatCapacity == true
        @parameters n kB θ
        @variables Tph(t)
        return nonlinear_phononheatcapacity(Tph,n,kB,θ)
    else
        return mp.Cph
    end
end

function nonlinear_phononheatcapacity(Tph::Float64,n::Float64,kB::Float64,θ::Float64)
    int(u,p) = nonlinear_phononheatcapacity_int(u)
    return 9*n*kB*(Tph/θ)^3*solve(IntegralProblem(int,0.0,θ/Tph),HCubatureJL(initdiv=2);abstol=1e-3,reltol=1e-3)
end
@register_symbolic nonlinear_phononheatcapacity(Tph::Num,n::Num,kB::Num,θ::Num)

function nonlinear_phononheatcapacity_int(u)
    return u^4*exp(u)/(exp(u)-1)^2
end

sim = define_simulation_settings(nlchempot=true,nlelecphon=true,nlelecheat=true)
mp = define_material_parameters(extcof=1.27e-8,gamma=71,debye=165,noatoms=59,plasma=2.1357,thermalcond=320.0,
elecperatom=1,eleceffmass=1.1,dos="DOS/Au_DOS.dat",secmomspecfun=29e6,elecphon=2.3e16,ballistic=0.0,cph=0.0)
@named Phonon_temp=t_phonon_factory(mp,sim)
function t_phonon_factory(mp::MaterialParameters,sim::SimulationSettings;name)
    @variables Tph(t) Tel(t)
    @named dTph = t_phonon_template()
    connections=[dTph.Source ~ t_phonon_sourceterm(sim),
                 dTph.HeatCapacity ~ t_phonon_heatcapacity(mp,sim),
                 dTph.ElecPhon ~ -1*t_electron_phononcoupling(mp,sim),
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
    return 9*n*kB*(Tph/θ)^3*solve(IntegralProblem(int,0.0,θ/Tph),HCubatureJL(initdiv=2);abstol=1e-3,reltol=1e-3).u
end
@register_symbolic nonlinear_phononheatcapacity(Tph::Num,n::Num,kB::Num,θ::Num)

function nonlinear_phononheatcapacity_int(u)
    return u^4*exp(u)/(exp(u)-1)^2
end
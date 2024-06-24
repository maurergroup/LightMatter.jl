function phonontemperature_factory(sim::SimulationSettings)
    HeatCapacity = phonontemperature_heatcapacity(sim)
    ElecPhon = Expr(:call,:*,-1,electronphonon_coupling(sim))
    Source = phonontemperature_source(sim)
    return build_phonontemperature(Source,ElecPhon,HeatCapacity)
end

function build_phonontemperature(Source,ElecPhon,HeatCapacity)
    return Expr(:call,:/,Expr(:call,:+,Source,ElecPhon),HeatCapacity)
end

function phonontemperature_heatcapacity(sim::SimulationSettings)
    if sim.ParameterApprox.PhononHeatCapacity == true
        return :(nonlinear_phononheatcapacity(Tph,mp.n,cons.kB,mp.θ))
    else
        return :(mp.Cph)
    end
end

function nonlinear_phononheatcapacity(Tph::Real,n::Real,kB::Real,θ::Real)
    int(u,p) = nonlinear_phononheatcapacity_int(u)
    return 9*n*kB*(Tph/θ)^3*solve(IntegralProblem(int,0.0,θ/Tph),HCubatureJL(initdiv=2);abstol=1e-3,reltol=1e-3).u
end
"""
    The integrand for the non-linear phononic heat capacity. Currently it is out-of-place.
"""
function nonlinear_phononheatcapacity_int(u::Real)
    return u^4*exp(u)/(exp(u)-1)^2
end

function phonontemperature_source(sim::SimulationSettings)
    if sim.Systems.NonEqElectrons == true
        if sim.Interactions.ElectronPhonon == true
            return :(electronphononscattering())
        else
            return 0.0
        end
    else
        return 0.0
    end
end
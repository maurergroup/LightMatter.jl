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
    int(u,p) = u^4*exp(u)/(exp(u)-1)^2
    prob = IntegralProblem(int,(0.0,θ/Tph))
    return 9*n*kB*(Tph/θ)^3*solve(prob,HCubatureJL(initdiv=2);abstol=1e-3,reltol=1e-3).u
end

function phonontemperature_source(sim::SimulationSettings)
    if sim.Systems.NonEqElectrons == true
        if sim.Interactions.ElectronPhonon == true
            return :(neqelectron_phonontransfer(fneq,mp.egrid,mp.τep,mp.u0,mp.FE,mp.DOS))
        else
            return 0.0
        end
    else
        return 0.0
    end
end

function neqelectron_phonontransfer(fneq,egrid,τep,u0,FE,DOS)
    spl=get_interpolate(egrid,fneq./τep)
    return get_internalenergyspl(spl,DOS,u0,FE)
end
function electrontemperature_factory(sim::SimulationSettings,laser::Expr,dim::Dimension)
    if sim.Systems.NonEqElectrons == false
        HeatCapacity = electrontemperature_heatcapacity(sim)
        ElecPhon = electronphonon_coupling(sim)
        Source = laser
        Spatial = electrontemperature_conductivity(dim::Dimension)
        return build_electronTTM(Source,Spatial,ElecPhon,HeatCapacity)
    elseif sim.Systems.NonEqElectrons == true
        Δu = :(elec_energychange(mp.egrid,relax_dis,μ,mp.DOS,u0))
        return build_athemelectron(Δu) 
    end
end

function build_electronTTM(Source,Spatial,ElecPhon,HeatCapacity)
    return Expr(:call,:/,Expr(:call,:+,Source,Spatial,ElecPhon),HeatCapacity)
end

function electrontemperature_heatcapacity(sim::SimulationSettings)
    if sim.ParameterApprox.ElectronHeatCapacity == true
        return :(nonlinear_electronheatcapacity(cons.kB,Tel,μ,mp.DOS))
    else
        return :(mp.γ*Tel)
    end
end

function nonlinear_electronheatcapacity(kB::Real,Tel::Real,μ::Real,DOS::Interpolations.Extrapolation)
    int(u,p) = dFDdT(kB,Tel,μ,u)*DOS(u)*u
    prob = IntegralProblem(int,(μ-(60*Tel/10000),μ+(60*Tel/10000)))
    return solve(prob,HCubatureJL(initdiv=5);reltol=1e-3,abstol=1e-3).u
end

function electronphonon_coupling(sim)
    if sim.Interactions.ElectronPhonon == true
        if sim.ParameterApprox.ElectronPhononCoupling==true
            return :(nonlinear_electronphononcoupling(cons.hbar,cons.kB,mp.λ,mp.DOS,Tel,μ,Tph,0.0))
        else
            return :(-mp.g*(Tel-Tph))
        end
    else
        return 0.0
    end
end

function nonlinear_electronphononcoupling(hbar::Real,kB::Real,λ::Real,DOS::Interpolations.Extrapolation,Tel::Real,μ::Real,Tph::Real,FE::Real)
    prefac=pi*kB*λ/DOS(μ)/hbar
    int(u,p) = DOS(u)^2*-dFDdE(kB,Tel,μ,u)
    prob = IntegralProblem(int,(μ-(60*Tel/10000),μ+(60*Tel/10000)))
    g=prefac.*solve(prob,HCubatureJL(initdiv=5);reltol=1e-3,abstol=1e-3).u
    return -g*(Tel-Tph)
end

function electrontemperature_conductivity(dim::Dimension)
    if typeof(dim) == Homogenous
        return 0.0
    else 
        return :(electronicthermalconductivity())
    end
end

function build_athemelectron(Δu)
    return :( 1/(c_T(μ,Tel,mp.DOS,cons.kB)*p_μ(μ,Tel,mp.DOS,cons.kB)-p_T(μ,Tel,mp.DOS,cons.kB)*c_μ(μ,Tel,mp.DOS,cons.kB))*(p_μ(μ,Tel,mp.DOS,cons.kB)*$Δu-c_μ(μ,Tel,mp.DOS,cons.kB)*Δn))
end

function elec_energychange(egrid,relax_dis,μ,DOS,u0)
    spl = get_interpolate(egrid,relax_dis)
    return get_internalenergyspl(μ,spl,DOS,u0)
end
function electrontemperature_factory(sim::SimulationSettings,laser::Expr,dim::Dimension)
    HeatCapacity = electrontemperature_heatcapacity(sim)
    ElecPhon = electronphonon_coupling(sim)
    Source = electrontemperature_source(sim,laser)
    Spatial = electrontemperature_conductivity(dim::Dimension)
    return build_electrontemperature(Source,Spatial,ElecPhon,HeatCapacity)
end

function build_electrontemperature(Source,Spatial,ElecPhon,HeatCapacity)
    return Expr(:call,:/,Expr(:call,:+,Source,Spatial,ElecPhon),HeatCapacity)
end

function electrontemperature_heatcapacity(sim::SimulationSettings)
    if sim.ParameterApprox.ElectronHeatCapacity == true
        return :(nonlinear_electronheatcapacity(cons.kB,Tel,μ-mp.FE,mp.DOS_nil))
    else
        return :(mp.γ*Tel)
    end
end

function nonlinear_electronheatcapacity(kB::Real,Tel::Real,μ::Real,DOS::Spline1D)
    p=(kB,Tel,μ,DOS)
    int(u,p) = electronheatcapacity_int(u,p)
    return solve(IntegralProblem(int,(μ-(6*Tel/10000),μ+(6*Tel/10000)),p),HCubatureJL(initdiv=50);reltol=1e-3,abstol=1e-3).u
end
"""
    The integrand for the non-linear electronic heat capacity using a parameter tuple 
    p=(kB, Tel, μ, DOS). The integrand is currently out-of-place.
"""
function electronheatcapacity_int(u::Real,p::Tuple{Real,Real,Real,Spline1D})
    return dFDdT(p[1],p[2],p[3],u)*p[4](u)*u
end

function electronphonon_coupling(sim)
    if sim.Interactions.ElectronPhonon == true
        if sim.ParameterApprox.ElectronPhononCoupling==true
            return :(nonlinear_electronphononcoupling(cons.hbar,cons.kB,mp.λ,mp.DOS_nil,Tel,μ-mp.FE,Tph,0.0))
        else
            return :(-mp.g*(Tel-Tph))
        end
    else
        return 0.0
    end
end

function nonlinear_electronphononcoupling(hbar::Real,kB::Real,λ::Real,DOS::Spline1D,Tel::Real,μ::Real,Tph::Real,FE::Real)
    prefac=pi*kB*λ/DOS(FE)/hbar
    p=(kB,Tel,μ,DOS)
    int(u,p) = electronphononcoupling_int(u,p)
    g=prefac.*solve(IntegralProblem(int,(μ-(6*Tel/10000),μ+(6*Tel/10000)),p),HCubatureJL(initdiv=10);reltol=1e-3,abstol=1e-3).u
    return -g*(Tel-Tph)
end
"""
    The integrand for the non-constant electron-phonon coupling parameter using a parameter tuple 
    p=(kB, Tel, μ, DOS). The integrand is currently out-of-place.
"""
function electronphononcoupling_int(u::Real,p::Tuple{Real,Real,Real,Spline1D})
    return p[4](u)^2*-dFDdE(p[1],p[2],p[3],u)
end

function electrontemperature_source(sim::SimulationSettings,laser::Expr)
    if sim.Systems.NonEqElectrons == true
        if sim.Interactions.ElectronElectron == true
            return neqelectron_electrontransfer()
        else
            return 0.0
        end
    else
        return laser
    end
end

function electrontemperature_conductivity(dim::Dimension)
    if typeof(dim) == Homogenous
        return 0.0
    else 
        return :(electronicthermalconductivity())
    end
end
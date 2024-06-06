
abstract type ElectronicTemperature <: Templates end

function (::ElectronicTemperature)(;name)
    @variables Tel(t) Source(t) ElecPhon(t) HeatCapacity(t) Spatial(t)

    eqs = D(Tel) ~ (Source .+ Spatial .+ ElecPhon)./HeatCapacity

    ODESystem(eqs,t;name)
end

function T_Electron_Factory(sim::SimulationSettings,laser::Num)
    @named dTel = (::ElectronicTemperature)
    source = T_Electron_Source(sim,laser)
    coupling = T_Electron_Coupling(sim)
    spatial = T_Electron_Transport(sim,grid)
    heatcapacity = T_Electron_HeatCapacity()
end

function T_Electron_Source(sim::SimulationSettings,laser::Num,mp::MaterialParameters,
    misc::Simulation_Misc)
    if sim.Interactions.Nonequilibrium_electrons == true
        return neq_elec_int_en(neqFD,misc,μ,Tel,gc,mp)
    elseif sim.Interactions.Nonequilibrium_electrons == true
        return laser
    end
end

function T_Electron_Coupling(sim::SimulationSettings)
    if sim.Interactions.ElectronPhonon == true
        return -1*ElectronPhononCoupling()
    elseif sim.Interactions.ElectronPhonon ==false
        return 0.0
    end
end

function ElectronPhononCoupling()
    @variables g Tel Tph
    return g*(Tel-Tph)
end

function T_Electron_HeatCapacity(sim::SimulationSettings)
    if sim.Nonlinear.ElectronHeatCapacity == true
        return Nonlinear_ElectronHeatCapacity()
    else
        return Linear_ElectronHeatCapacity()
    end
end

function Nonlinear_ElectronHeatCapacity()
    int(u,p) = DOS(u) * u
end

function Linear_ElectronHeatCapacity()
end

function T_Electron_Transport(sim::SimulationSettings,grid::Dimension)
    if typeof(grid) == Homogenous
        return 0.0
    elseif grid.Dim != 0
        return Thermal_Electron_Transport()
    end
end

function FermiDirac(E::Float64,μ::Float64,Tel,gc)::Float64
    return 1/(exp((E-μ)/(gc.kB*Tel))+1)
end




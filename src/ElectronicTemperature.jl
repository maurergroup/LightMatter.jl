"""
    electrontemperature_factory(sim::Simulation, laser::Expr)
    
    Assembles expresssion for how the electronic temperature should be propgated through time.

    # Arguments
    - 'sim': Simulation settings and parameters
    - 'laser': Expression for laser energy as a funciton of time

    # Returns
    - Expression for the time-propagation of the electronic temperature subsystem
"""
function electrontemperature_factory(sim::Simulation, laser::Expr)
    if sim.electronictemperature.AthermalElectron_ElectronCoupling == false
        HeatCapacity = electrontemperature_heatcapacity(sim)
        ElecPhon = electronphonon_coupling(sim)
        return build_electronTTM(sim, laser, ElecPhon, HeatCapacity)
    elseif sim.electronictemperature.AthermalElectron_ElectronCoupling == true
        Δu = athem_electempenergychange(sim)
        return build_athemelectron(Δu) 
    end
end
"""
    build_electronTTM(sim::Simulation, Source::Expr, ElecPhon::Expr, HeatCapacity::Expr)
    
    Builds the differential equation (expression) for a Two-Temperature like system.

    # Arguments
    - 'sim': Simulation settings and parameters
    - 'Source': Expression for the incoming energy source (typically a laser)
    - 'ElecPhon': Expression for the electron-phonon coupling
    - 'HeatCapacity': Expression for calculating the heat capacity

    # Returns
    - Expression for the time evolution of a two-temperature model electronic temperature
"""
function build_electronTTM(sim::Simulation, Source::Expr, ElecPhon::Expr, HeatCapacity::Expr)
    args = Union{Expr, Symbol, Float64}[Source, ElecPhon]
    if sim.electronictemperature.Conductivity == true
        push!(args, :Tel_cond)
    end
    return Expr(:call, :/, Expr(:call, :+, args...), HeatCapacity)
end
"""
    electrontemperature_heatcapacity(sim::Simulation)
    
    Determines the expression for the electronic temperature's heat capacity.
    Currently implemented:
    - :linear : Specific heat of electrons(γ) multiplied by the current temperature
    - :nonlinear : Calculated from the density-of-states of the system

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Expression for the electronic temperature's heat capacity
"""
function electrontemperature_heatcapacity(sim::Simulation)
    if sim.electronictemperature.ElectronicHeatCapacity == :nonlinear
        return :(LightMatter.nonlinear_electronheatcapacity(Tel, μ, DOS))
    elseif sim.electronictemperature.ElectronicHeatCapacity == :linear
        return :(sim.electronictemperature.γ * Tel)
    elseif sim.electronictemperature.ElectronicHeatCapacity == :constant
        return :(sim.electronictemperature.γ)
    end
end
"""
    nonlinear_electronheatcapacity(Tel::Float64, μ::Float64, DOS::spl, egrid::Vector{Float64})
    
    Calculates non-linear electronic bath heat capacity. A more accurate method than the
    linear form.

    # Arguments
    - 'Tel': Temperature of the electronic bath
    - 'μ': Chemical potential of the electronic bath
    - 'DOS': Density-of-states of the system
    - 'egrid': Energy grid the simulation is solved over

    # Returns
    - The current heat capacity of the electronic thermal bath
"""
function nonlinear_electronheatcapacity(Tel, μ, DOS)
    int(u,p) = dFDdT(Tel, μ, u) * DOS(u) * u
    prob = IntegralProblem(int, -8*Constants.kB*Tel, 8*Constants.kB*Tel)
    return solve(prob, HCubatureJL(initdiv=25, buffer=true), abstol=1e-5, reltol=1e-5).u
end
"""
    electronphonon_coupling(sim::Simulation)
    
    Determines the expression for the coupling between an electronic and phononic thermal bath.
    Currently implemented:
    - :variable : Calculates the electron-phonon coupling parameter from the density-of-states of the system
    - :constant : Uses a constant value for the electron-phonon coupling parameter (g)

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Expression for the electron-phonon thermal bath energy transfer
"""
function electronphonon_coupling(sim::Simulation)
    if sim.electronictemperature.Electron_PhononCoupling == true
        if sim.electronictemperature.ElectronPhononCouplingValue == :variable
            return :(LightMatter.variable_electronphononcoupling(sim.electronictemperature.λ, sim.electronictemperature.ω, DOS, Tel, μ, Tph))
        elseif sim.electronictemperature.ElectronPhononCouplingValue == :constant
            return :(-sim.electronictemperature.g*(Tel-Tph))
        end
    else
        return 0.0
    end
end
"""
    variable_electronphononcoupling(λ::Float64, ω::Float64, DOS::spl, Tel::Float64, μ::Float64, Tph::Float64, egrid::Vector{Float64})
    
    Calculates the non-linear electron phonon coupling parameter and subsequent energy flow from the density-of-states
    of the system. More accurate than a constant value. 
    The expression can be found in Z. Lin, L. V. Zhigilei and V. Celli, Phys. Rev. B, 2008, 77, 075133.

    # Arguments
    - 'λ': Electron-phonon mass enhancement parameter
    - 'ω': Second moment of the phonon spectrum
    - 'DOS': Density-of-states of the system
    - 'Tel': Temperature of the electronic bath
    - 'μ': Chemical potential of the electronic bath
    - 'Tph': Temperature of the phonon bath
    - 'egrid': Energy grid the simulation is solved over

    # Returns
    - Energy flow between an electronic and phononic bath with a calculate g parameter
"""
function variable_electronphononcoupling(λ, ω, DOS, Tel, μ, Tph)
    prefac=pi * Constants.kB * λ * ω / DOS(μ) / Constants.ħ
    int(u,p) = DOS(u)^2 *-dFDdE(Tel, μ, u) * prefac
    prob = IntegralProblem(int, -8*Constants.kB*Tel, 8*Constants.kB*Tel)
    sol = solve(prob, HCubatureJL(initdiv=25, buffer=true), abstol=1e-5, reltol=1e-5).u
    return -sol * (Tel-Tph)
end
"""
    build_athemelectron(Δu::Expr)
    
    Builds the differential equation (expression) for an AthEM thermal bath. The difference to the TTM
    is due to accounting for the changing particle number

    # Arguments
    - 'Δu': Expression for the change in internal energy of the electronic bath

    # Returns
    - Expression for the time evolution of a, AthEM electronic temperature
"""
function build_athemelectron(Δu::Expr)
    return :( 1 / (LightMatter.c_T(μ,Tel,DOS,sim.structure.egrid)*LightMatter.p_μ(μ,Tel,DOS,sim.structure.egrid)
    - LightMatter.p_T(μ,Tel,DOS,sim.structure.egrid)*LightMatter.c_μ(μ,Tel,DOS,sim.structure.egrid)) *
    (LightMatter.p_μ(μ,Tel,DOS,sim.structure.egrid)*$Δu - LightMatter.c_μ(μ,Tel,DOS,sim.structure.egrid)*du.noe[i]))
end
"""
    athem_electempenergychange(sim::Simulation)
    
    Builds an expression for the change in internal energy of an AthEM electronic bath

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Expression for the change of internal energy of an AthEM electronic bath
"""
function athem_electempenergychange(sim::Simulation)
    args = Vector{Union{Expr,Symbol}}(undef, 0)
    push!(args, :(LightMatter.elec_energychange(sim.structure.egrid, relax_dis, DOS)))
    if sim.phononictemperature.Enabled == true
       push!(args, electronphonon_coupling(sim::Simulation))
    end
    if sim.electronictemperature.Conductivity == true
        push!(args, :(Tel_cond))
    end
    return Expr(:call, :+, args...)
end
"""
    elec_energychange(egrid::Vector{Float64}, relax_dis::Vector{Float64}, DOS::spl)
    
    Calculates the energy change of the thermal bath due to non-equilibrium-equilibrium
    electron scattering.

    # Arguments
    - 'egrid: Energy grid the simulation is solved over
    - 'relax_dis': The change in distribution due to the electron-electron scattering
    - 'DOS': Density-of-states of the system

    # Returns
    - The change in the internal energy of the thermal system due to e-e scattering
"""
function elec_energychange(egrid, relax_dis, DOS)
    return LightMatter.get_internalenergy(relax_dis, DOS, egrid)
end
"""
    electrontemperature_conductivity!(Tel::Vector{Float64}, κ::Union{Float64,Vector{Float64}}, dz::Float64, Tph::Vector{Float64}, cond::Vector{Float64})
    
    Calculates the change in temperature due to diffusive transport of energy through the system. Uses the boundary conditions of
    setting dTel / dz to 0.0.

    # Arguments
    - 'Tel': Temperature of the electronic bath
    - 'κ': Electronic thermal conductivity at room temperature
    - 'dz': z-grid spacing
    - 'Tph': Temperature of the phonon bath
    - 'cond': Vector to store the change in electronic temperature

    # Returns
    - Updates the cond vector with the change in electronic temperature at each grid point
"""
function electrontemperature_conductivity!(Tel, κ, dz, Tph, cond)
    cond = get_tmp(cond, Tel[1])
    K = κ.*Tel./Tph

    for i in 2:length(K)-1
        K_plus = 1 / 2 * (K[i+1] + K[i])
        K_minus = 1 / 2 * (K[i] + K[i-1])
        cond[i] = (K_plus*(Tel[i+1] - Tel[i]) - K_minus*(Tel[i] - Tel[i-1])) / dz^2
    end

    K_plus1 = 1 / 2 * (K[2] + K[1])
    cond[1] = (K_plus1*(Tel[2] - Tel[1]) ) / dz^2
    K_plusend = 1 / 2 * (K[end] + K[end-1])
    cond[end] = -(K_plusend*(Tel[end] - Tel[end-1]) ) / dz^2
end
"""
    depthderivative!(vec::Vector{Float64}, dz::Float64, Diff::Vector{Float64})
    
    Calculates a 2nd order finite difference method of a vector along a grid with spacing dz.
    Uses central difference in the middle and forward(reverse) difference at the top(end) of the vector

    # Arguments
    - 'vec': Vector the finite difference is being performed over
    - 'dz': z-grid spacing
    - 'Diff': Vector to store the finite difference result

    # Returns
    - Updates the diff vector with the finite difference result
"""
function depthderivative!(vec, dz, Diff)
    for i in 2:length(vec)-1
        Diff[i]=(vec[i+1]-vec[i-1])/(2*dz)
    end
    Diff[1] = (vec[2]-vec[1])/dz
    Diff[end] = (vec[end]-vec[end-1])/dz
end
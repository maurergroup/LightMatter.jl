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
    args = Union{Expr, Symbol, Number}[Source, ElecPhon]
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
        return :(Lightmatter.nonlinear_electronheatcapacity(Tel, μ, DOS, sim.structure.egrid))
    elseif sim.electronictemperature.ElectronicHeatCapacity == :linear
        return :(sim.electronictemperature.γ * Tel)
    end
end
"""
    nonlinear_electronheatcapacity(Tel::Number, μ::Number, DOS::spl, egrid::Vector{<:Number})
    
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
function nonlinear_electronheatcapacity(Tel, μ, DOS, egrid)
    return extended_Bode(dFDdT(Tel,μ,egrid).*DOS(egrid).*egrid, egrid)
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
            return :(Lightmatter.nonlinear_electronphononcoupling(sim.electronictemperature.λ, sim.electronictemperature.ω, DOS, Tel, μ, Tph, sim.structure.egrid))
        elseif sim.electronictemperature.ElectronPhononCouplingValue == :constant
            return :(-sim.electronictemperature.g*(Tel-Tph))
        end
    else
        return 0.0
    end
end
"""
    nonlinear_electronphononcoupling(λ::Number, ω::Number, DOS::spl, Tel::Number, μ::Number, Tph::Number, egrid::Vector{<:Number})
    
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
    - Energy flow between an electornic and phononic bath with a calculate g parameter
"""
function nonlinear_electronphononcoupling(λ, ω, DOS, Tel, μ, Tph, egrid)
    prefac=pi * Constants.kB * λ * ω / DOS(μ) / Constants.ħ
    g=prefac .* extended_Bode(DOS(egrid).^2 .*-dFDdE(Tel, μ, egrid), egrid)
    return -g * (Tel-Tph)
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
    return :( 1 / (Lightmatter.c_T(μ,Tel,DOS,sim.structure.egrid)*Lightmatter.p_μ(μ,Tel,DOS,sim.structure.egrid)
    - Lightmatter.p_T(μ,Tel,DOS,sim.structure.egrid)*Lightmatter.c_μ(μ,Tel,DOS,sim.structure.egrid)) *
    (Lightmatter.p_μ(μ,Tel,DOS,sim.structure.egrid)*$Δu - Lightmatter.c_μ(μ,Tel,DOS,sim.structure.egrid)*Δn))
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
    push!(args, :(Lightmatter.elec_energychange(sim.structure.egrid, relax_dis, DOS)))
    if sim.phononictemperature.Enabled == true
       push!(args, electronphonon_coupling(sim::Simulation))
    end
    if sim.electronictemperature.Conductivity == true
        push!(args, :(Tel_cond))
    end
    return Expr(:call, :+, args...)
end
"""
    elec_energychange(egrid::Vector{<:Number}, relax_dis::Vector{<:Number}, DOS::spl)
    
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
    return Lightmatter.get_internalenergy(relax_dis, DOS, egrid)
end
"""
    electrontemperature_conductivity!(Tel::Vector{<:Number}, κ::Union{Number,Vector{<:Number}}, dz::Number, Tph::Vector{<:Number}, cond::Vector{<:Number})
    
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
    depthderivative!(vec::Vector{<:Number}, dz::Number, Diff::Vector{<:Number})
    
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

function harmonic_average_K(K)
    new_K = zeros(length(K))
    new_K[1:end-1] = 2 .* K[2:end] .* K[1:end-1] ./ (K[2:end] .+ K[1:end-1])
    new_K[end] = 2 * K[end] * K[end-1] ./ (K[end-1] + K[end])
    return new_K
end
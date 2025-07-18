"""
    phonontemperature_factory(sim::Simulation)
    
    Assembles expresssion for how the phononic temperature should be propgated through time.

    # Arguments
    - 'sim': Simulation settings and parameters
    # Returns
    - Expression for the time-propagation of the phononic temperature subsystem
"""
function phonontemperature_factory(sim::Simulation)
    HeatCapacity = phonontemperature_heatcapacity(sim)
    ElecPhon = Expr(:call, :*, -1, electronphonon_coupling(sim))
    Source = phonontemperature_source(sim)
    return build_phonontemperature(sim, Source, ElecPhon, HeatCapacity)
end
"""
    build_phonontemperature(sim::Simulation, Source::Expr, ElecPhon::Expr, HeatCapacity::Expr)
    
    Builds the differential equation (expression) for the phonon bath

    # Arguments
    - 'sim': Simulation settings and parameters
    - 'Source': Expression for energy input from athermal electrons
    - 'ElecPhon': Expression for the thermal electron-phonon coupling
    - 'HeatCapacity': Expression for calculating the heat capacity

    # Returns
    - Expression for the time evolution of a phnonic thermal bath
"""
function build_phonontemperature(sim::Simulation, Source::Union{Expr,Float64}, ElecPhon::Expr, HeatCapacity::Expr)
    args = Union{Expr,Symbol,Float64}[Source,ElecPhon]
    if sim.phononictemperature.Conductivity == true
        push!(args,:Tph_cond)
    end
    return Expr(:call,:./,Expr(:call,:+,args...),HeatCapacity)
end
"""
    phonontemperature_heatcapacity(sim::Simulation)
    
    Determines the expression for the phononic temperature's heat capacity.
    Currently implemented:
    - :constant : Constant phonon heat capacity (Cph)
    - :nonlinear : Calculated from Simpson's rule

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Expression for the phononic temperature's heat capacity
"""
function phonontemperature_heatcapacity(sim::Simulation)
    if sim.phononictemperature.PhononicHeatCapacity == :nonlinear
        return :(LightMatter.nonlinear_phononheatcapacity(Tph, sim.phononictemperature.n, sim.phononictemperature.θ))
    elseif sim.phononictemperature.PhononicHeatCapacity == :constant
        return :(sim.phononictemperature.Cph)
    end
end
"""
    nonlinear_phononheatcapacity(Tph::Float64, n::Float64, θ::Float64)
    
    Calculates non-linear phononic bath heat capacity. A more accurate method than the
    constant form.

    # Arguments
    - 'Tph': Temperature of the phononic bath
    - 'n': Float64 of atoms per nm³
    - 'θ': Debye temperature of the system

    # Returns
    - The current heat capacity of the phononic thermal bath
"""
function nonlinear_phononheatcapacity(Tph, n, θ)
    int(u,p) = u^4 * exp(u) / (exp(u)-1)^2
    prob = IntegralProblem(int, (0.0, θ/Tph))
    return 9*n*Constants.kB*(Tph/θ)^3 * solve(prob, HCubatureJL(initdiv=10); abstol=1e-5, reltol=1e-5).u
end

#= function nonlinear_phononheatcapacity(Tph::ForwardDiff.Dual, n, θ)
    temp = ForwardDiff.value(Tph)
    int(u,p) = u^4 * exp(u) / (exp(u)-1)^2
    prob = IntegralProblem(int, (0.0, θ/temp))
    return 9*n*Constants.kB*(Tph/θ)^3 * solve(prob, HCubatureJL(initdiv=10); abstol=1e-5, reltol=1e-5).u
end =#
"""
    phonontemperature_source(sim::Simulation)
    
    Determines the expresision of any additional source terms into the phonon bath. Currently this
    is just athermal electron-phonon coupling
    Currently implemented:

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Expression for any source term into the phonon bath
"""
function phonontemperature_source(sim::Simulation)
    if sim.phononictemperature.AthermalElectron_PhononCoupling == true
        τep = phonon_relaxationtime(sim::Simulation)
        return :(LightMatter.neqelectron_phonontransfer(fneq, sim.structure.egrid, $τep, DOS))
    else
        return 0.0
    end
end
"""
    neqelectron_phonontransfer(fneq::Vector{Float64}, egrid::Vector{Float64}, τep::Float64, DOS::spl)
    
    Calculates energy input into the phonon bath due to non-equilibrium electron-phonon scattering

    # Arguments
    - 'fneq': Non-equilibrium electron distribution
    - 'egrid': Energy grid the distributions are solved on
    - 'τep': Lifetime of athermal electrons due to electron-phonon scattering
    - 'DOS': Density-of-states of the system

    # Returns
    - Value of the change in the phonon internal energy
"""
function neqelectron_phonontransfer(fneq, egrid, τep, DOS)
    return LightMatter.get_internalenergy(fneq./τep,DOS,egrid)
end
"""
    phonontemperature_conductivity!(Tph::Vector{Float64}, κ::Union{Float64,Vector{Float64}}, dz::Float64, cond::Vector{Float64})
    
    Calculates thermal phonon conductivity due to diffusive transport

    # Arguments
    - 'Tph': Phonon bath temperature
    - 'κ': Phonon thermal conductivity (assumed constant)
    - 'dz': Spacing of the spatial grid
    - 'cond': Vector to store the temperature change

    # Returns
    - Updates cond with the change in temperature at each z-grid point
"""
function phonontemperature_conductivity!(Tph, κ, dz, cond)
    depthderivative!(Tph, dz, cond)
    cond[1] = 0.0
    cond[end] = 0.0
    depthderivative!((cond.*κ), dz, cond)
end

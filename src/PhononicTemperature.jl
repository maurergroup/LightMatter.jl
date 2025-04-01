"""
    phonontemperature_factory(sim::Simulation)
    This function takes the Simulation struct and returns an assembled
    expression for how the phononic temperature should evolve over time. 
"""
function phonontemperature_factory(sim::Simulation)
    HeatCapacity = phonontemperature_heatcapacity(sim)
    ElecPhon = Expr(:call,:*,-1,electronphonon_coupling(sim))
    Source = phonontemperature_source(sim)
    return build_phonontemperature(sim,Source,ElecPhon,HeatCapacity)
end
"""
    build_phonontemperature(sim::Simulation,Source::Union{Expr,Real},ElecPhon::Expr,HeatCapacity::Expr)
    Builds a combine expression of the expressions for the energy input, electorn-phonon coupling
    and the phononic heat capacity. The first two terms are summed to find
    the change in internal energy of the thermal system and divided by the heat capacity to change the 
    internal energy into a temperature.
"""
function build_phonontemperature(sim::Simulation,Source::Union{Expr,Real},ElecPhon::Expr,HeatCapacity::Expr)
    args = Union{Expr,Symbol,Real}[Source,ElecPhon]
    if sim.phononictemperature.Conductivity == true
        push!(args,:Tph_cond)
    end
    return Expr(:call,:./,Expr(:call,:+,args...),HeatCapacity)
end
"""
    phonontemperature_heatcapacity(sim::Simulation)
    Returns the expression for how the phononic heat capacity should be 
    calculated. This can be done as a constant (mp.Cph) or via a nonlinear
    relationship. The keyword to be set in define_simulation_settings is nlphonheat = true.
"""
function phonontemperature_heatcapacity(sim::Simulation)
    if sim.phononictemperature.PhononicHeatCapacity == :nonlinear
        return :(Lightmatter.nonlinear_phononheatcapacity(Tph,sim.phononictemperature.n,sim.phononictemperature.θ))
    elseif sim.phononictemperature.PhononicHeatCapacity == :constant
        return :(sim.phononictemperature.Cph)
    end
end
"""
    nonlinear_phononheatcapacity(Tph::Real,n::Real,θ::Real)
    This is the function that returns the nonlinear variant of the phononic heat capacity.
    It is found using Simpson's rule.
"""
function nonlinear_phononheatcapacity(Tph::Real,n::Real,θ::Real)
    int(u,p) = u^4*exp(u)/(exp(u)-1)^2
    prob = IntegralProblem(int,(0.0,θ/Tph))
    return 9*n*Constants.kB*(Tph/θ)^3*solve(prob,HCubatureJL(initdiv=10);abstol=1e-5,reltol=1e-5).u
end
"""
    phonontemperature_source(sim::Simulation)
    Selects the energy source into the phonons. This does not include thermal electron-
    phonon coupling which is defined seperately but instead returns the input from the athermal electrons
    in AthEM to the phonon system for example. 
"""
function phonontemperature_source(sim::Simulation)
    if sim.phononictemperature.AthermalElectron_PhononCoupling == true
        τep = phonon_relaxationtime(sim::Simulation)
        return :(Lightmatter.neqelectron_phonontransfer(fneq,sim.structure.egrid,$τep,DOS))
    else
        return 0.0
    end
end
"""
    neqelectron_phonontransfer(fneq::Vector{<:Real},egrid::Vector{<:Real},τep::Real,DOS::spl)
    Calculates the energy transfer from the athermal electrons in AthEM to the thermal phonons.
"""
function neqelectron_phonontransfer(fneq::Vector{<:Real},egrid::Vector{<:Real},τep::Real,DOS::spl)
    return Lightmatter.get_internalenergy(fneq./τep,DOS,egrid)
end
"""
    phonontemperature_conductivity!(Tph::Vector{<:Real},κ::Real,dz::Real,cond::Vector{<:Real})
    Calculates the thermal energy passing further into a bulk slab due to thermal conductivity of the phononic bath. 
    If the type of the Dimension struct is Homogeneous then there should be no conductivty and returns 0.0 at every time step. 
    The derivative of the temperature with respect to distance is set to 0 at the boundaries.
"""
function phonontemperature_conductivity!(Tph::Vector{<:Real},κ::Union{Real,Vector{<:Real}},dz::Real,cond::Vector{<:Real})
    depthderivative!(Tph,dz,cond)
    cond[1] = 0.0
    cond[end] = 0.0
    depthderivative!((cond.*κ),dz,cond)
end

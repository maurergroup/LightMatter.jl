"""
    phonontemperature_factory(sim::SimulationSettings)
    This function takes the SimulationSettings struct and returns an assembled
    expression for how the phononic temperature should evolve over time. 
"""
function phonontemperature_factory(sim::SimulationSettings)
    HeatCapacity = phonontemperature_heatcapacity(sim)
    ElecPhon = Expr(:call,:*,-1,electronphonon_coupling(sim))
    Source = phonontemperature_source(sim)
    return build_phonontemperature(Source,ElecPhon,HeatCapacity)
end
"""
    build_phonontemperature(Source::Expr,ElecPhon::Expr,HeatCapacity::Expr)
    Builds a combine expression of the expressions for the energy input, electorn-phonon coupling
    and the phononic heat capacity. The first two terms are summed to find
    the change in internal energy of the thermal system and divided by the heat capacity to change the 
    internal energy into a temperature.
"""
function build_phonontemperature(Source::Union{Expr,Real},ElecPhon::Expr,HeatCapacity::Expr)
    return Expr(:call,:./,Expr(:call,:+,Source,ElecPhon),HeatCapacity)
end
"""
    phonontemperature_heatcapacity(sim::SimulationSettings)
    Returns the expression for how the phononic heat capacity should be 
    calculated. This can be done as a constant (mp.Cph) or via a nonlinear
    relationship. The keyword to be set in define_simulation_settings is nlphonheat = true.
"""
function phonontemperature_heatcapacity(sim::SimulationSettings)
    if sim.ParameterApprox.PhononHeatCapacity == true
        return :(Lightmatter.nonlinear_phononheatcapacity(Tph,mp.n,cons.kB,mp.θ))
    else
        return :(mp.Cph)
    end
end
"""
    nonlinear_phononheatcapacity(Tph::Real,n::Real,kB::Real,θ::Real)
    This is the function that returns the nonlinear variant of the phononic heat capacity.
    It is found using Simpson's rule.
"""
function nonlinear_phononheatcapacity(Tph::Real,n::Real,kB::Real,θ::Real)
    int(u,p) = u^4*exp(u)/(exp(u)-1)^2
    prob = IntegralProblem(int,(0.0,θ/Tph))
    return 9*n*kB*(Tph/θ)^3*solve(prob,HCubatureJL(initdiv=2);abstol=1e-3,reltol=1e-3).u
end
"""
    phonontemperature_source(sim::SimulationSettings)
    Selects the energy source into the phonons. This does not include thermal electron-
    phonon coupling which is defined seperately but instead returns the input from the athermal electrons
    in AthEM to the phonon system for example. 
"""
function phonontemperature_source(sim::SimulationSettings)
    if sim.Systems.NonEqElectrons == true && sim.Interactions.ElectronPhonon == true
        return :(Lightmatter.neqelectron_phonontransfer(fneq,mp.egrid,mp.τep,DOS))
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
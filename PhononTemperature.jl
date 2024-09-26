"""
    This the base factory function that constructs the phononic temperature ODE. All variables and
    functionality for the phononic temperature ODE should be set up within this function call. For
    the electron-phonon coupling it currently inherits the same equation from that of the electronic 
    temperature. To-Do: Check performance if I just make a connection between the two electron-phonon
    coupling variables.
"""
function t_phonon_factory(mp::MaterialParameters,sim::SimulationSettings;name)
    @variables Tph(t) Tel(t)
    @named dTph = t_phonon_template()
    connections=[dTph.Source ~ t_phonon_sourceterm(sim),
                 dTph.HeatCapacity ~ t_phonon_heatcapacity(sim),
                 dTph.ElecPhon ~ -t_electron_phononcoupling(sim),
                 dTph.Tph ~ Tph]
    compose(ODESystem(connections,t;name),dTph)
end
"""
    This is the template for the phononic temperature ODE. It contains a source term where 
    energy injection functions are defined, a spatial term for thermal conducitivity within the
    phononic system. An electron-phonon coupling term for how the phononic temperature interacts
    with the electronic temperature and then the heat capcaity of the phononic system. All variables
    are then later replaced by their relative function for the current simulation. Any functionality
    that is unwanted during a simulation must be set to 0.0 during setup, not ignored.
"""
function t_phonon_template(;name)
    @variables Tph(t) Source(t) ElecPhon(t) HeatCapacity(t) Spatial(t)

    eqs = D(Tph) ~ (Source .+ ElecPhon)./HeatCapacity

    ODESystem(eqs,t;name)
end
"""
    Defines and returns the requested equation for how energy is inputed into the system.
    Currently the only options are either AthEM or the laser directly so checks whether non-
    equilibrium electrons are enabled and defaults to the AthEM method otherwise sets to 0.0.
"""
function t_phonon_sourceterm(sim::SimulationSettings)
    if sim.Systems.NonEqElectrons == true
        if sim.Interactions.ElectronPhonon == true
            @variables uep(t)
            return uep
        else
            return 0.0
        end
    else
        return 0.0
    end
end
"""
    Defines and returns the requested equation for the electronic heat capacity. Uses the Boolean
    flag to determine whether a constant parameter is used or Simpson's rule is used to evaluate 
    the heat capacity.
"""
function t_phonon_heatcapacity(sim::SimulationSettings)
    if sim.ParameterApprox.PhononHeatCapacity == true
        @parameters n kB θ
        @variables Tph(t)
        return nonlinear_phononheatcapacity(Tph,n,kB,θ)
    else
        @parameters Cph
        return Cph
    end
end
"""
    The non-linear phononic heat capacity equation using SImpson's Rule. It evaluates the 
    integrand using a serial solver.
"""
function nonlinear_phononheatcapacity(Tph::Real,n::Real,kB::Real,θ::Real)
    int(u,p) = nonlinear_phononheatcapacity_int(u)
    return 9*n*kB*(Tph/θ)^3*solve(IntegralProblem(int,0.0,θ/Tph),HCubatureJL(initdiv=2);abstol=1e-3,reltol=1e-3).u
end
@register_symbolic nonlinear_phononheatcapacity(Tph::Num,n::Num,kB::Num,θ::Num)
"""
    The integrand for the non-linear phononic heat capacity. Currently it is out-of-place.
"""
function nonlinear_phononheatcapacity_int(u::Real)
    return u^4*exp(u)/(exp(u)-1)^2
end
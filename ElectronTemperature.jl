"""
    This the base factory function that constructs the electronic temperature ODE. All variables and
    functionality for the electronic temperature ODE should be set up within this function call.
"""
function t_electron_factory(mp::MaterialParameters,sim::SimulationSettings,laser::Num;name)
    @variables Tel(t)
    @named dTel = t_elec_template()
    connections =[dTel.HeatCapacity ~ t_electron_heatcapacity(mp,sim),
                  dTel.ElecPhon ~ t_electron_phononcoupling(mp,sim),
                  dTel.Spatial ~ 0.0,
                  dTel.Source ~ t_electron_sourceterm(sim,laser),
                  dTel.Tel ~ Tel]
    compose(ODESystem(connections,t;name),dTel)
end
"""
    This is the template for the electonic temperature ODE. It contains a source term where 
    energy injection functions are defined, a spatial term for thermal conductiivity within the
    electronic system. An electron-phonon coupling term for how the electronic temperature interacts
    with the lattice temperature and then the heat capcaity of the electronic system. All variables
    are then later replaced by their relative function for the current simulation. Any functionality
    that is unwanted during a simulation must be set to 0.0 during setup, not ignored.
"""
function t_elec_template(;name)
    @variables Tel(t) Source(t) ElecPhon(t) HeatCapacity(t) Spatial(t)

    eqs = Tel ~ (Source .+ Spatial .+ ElecPhon)./HeatCapacity

    ODESystem(eqs,t;name)
end
"""
    Defines and returns the requested equation for the electronic heat capacity. Uses the Boolean
    flag to determine whether a linear approximation or non-linear term is utilised and sets up
    the relative equation.
"""
function t_electron_heatcapacity(mp::MaterialParameters,sim::SimulationSettings)
    if sim.ParameterApprox.ElectronHeatCapacity == true
        @parameters kB μ
        @variables  Tel(t)
        return nonlinear_electronheatcapacity(kB,Tel,μ,mp.DOS)
    else
        @parameters γ 
        @variables Tel(t)

        return γ*Tel
    end
end
"""
    The non-linear electronic heat capacity equation. It evaluates the integrand using a serial solver.
    The integrand limits are temperature dependent to account for errors due to discontinuties at high
    energy values. To-Do: Check limits work for many materials.
"""
function nonlinear_electronheatcapacity(kB,Tel,μ,DOS::Spline1D)
    p=(kB,Tel,μ,DOS)
    int(u,p) = electronheatcapacity_int(u,p)
    return solve(IntegralProblem(int,(μ-(6*Tel/10000),μ+(6*Tel/10000)),p),HCubatureJL(initdiv=10);reltol=1e-3,abstol=1e-3).u
end
@register_symbolic nonlinear_electronheatcapacity(kB::Num,Tel::Num,μ::Num,DOS::Spline1D)
"""
    The integrand for the non-linear electronic heat capacity using a parameter tuple 
    p=(kB, Tel, μ, DOS). The integrand is currently out-of-place.
"""
function electronheatcapacity_int(u::Real,p::Tuple{Real,Real,Real,Spline1D})
    return dFDdT(p[1],p[2],p[3],u)*p[4](u)*u
end
"""
    Defines and returns the requested equation for the electron-phonon coupling. Uses the Boolean
    flag to determine whether a constant or non-constant g value is used.
"""
function t_electron_phononcoupling(mp::MaterialParameters,sim::SimulationSettings)
    if sim.Interactions.ElectronPhonon == true
        if sim.ParameterApprox.ElectronPhononCoupling==true
            @parameters hbar kB λ μ
            @variables Tel(t) Tph(t)

            return nonlinear_electronphononcoupling(hbar,kB,λ,mp.DOS,Tel,μ,Tph)
        else
            @parameters g 
            @variables Tel(t) Tph(t)

            return -g*(Tel-Tph)
        end
    else
        return 0.0
    end
end
"""
    The non-constant electron-phonon coupling equation.. It evaluates the integrand using a serial 
    solver. The integrand limits are temperature dependent to account for errors due to 
    discontinuties at high energy values. To-Do: Check limits work for many materials.
"""
function nonlinear_electronphononcoupling(hbar::Real,kB::Real,λ::Real,DOS::Spline1D,Tel::Real,μ::Real,Tph::Real)
    prefac=pi*kB*λ/DOS(μ)/hbar
    p=(kB,Tel,μ,DOS)
    int(u,p) = electronphononcoupling_int(u,p)
    g=prefac.*solve(IntegralProblem(int,(μ-(6*Tel/10000),μ+(6*Tel/10000)),p),HCubatureJL(initdiv=10);reltol=1e-3,abstol=1e-3).u
    return -g*(Tel-Tph)
end
@register_symbolic nonlinear_electronphononcoupling(hbar::Num,kB::Num,λ::Num,DOS::Spline1D,Tel::Num,μ::Num,Tph::Num)
"""
    The integrand for the non-constant electron-phonon coupling parameter using a parameter tuple 
    p=(kB, Tel, μ, DOS). The integrand is currently out-of-place.
"""
function electronphononcoupling_int(u::Real,p::Tuple{Real,Real,Real,Spline1D})
    return p[4](u)^2*-dFDdE(p[1],p[2],p[3],u)
end
"""
    Defines and returns the requested equation for how energy is inputed into the system.
    Currently the only options are either AthEM or the laser directly so checks whether non-
    equilibrium electrons are enabled and defaults to the AthEM method otherwise uses the laser.
    The laser is currently brought into the function because using a variable and connecting later
    leads to un-assignable parameters. To-Do: Review once Simulation Builder is complete
"""
function t_electron_sourceterm(sim::SimulationSettings,laser::Num)
    if sim.Systems.NonEqElectrons == true
        if sim.Interactions.ElectronElectron == true
            @variables uee(t)
            return uee
        else
            return 0.0
        end
    else
        return laser
    end
end


"""
    This the base factory function that constructs the electronic temperature ODE. All variables and
    functionality for the electronic temperature ODE should be set up within this function call.
"""
function t_electron_factory(mp::MaterialParameters,sim::SimulationSettings,laser::Num;name)
    @variables Tel(t)
    if sim.Systems.NonEqElectrons==false
        @named dTel = ttm_elec_template()
        connections =[dTel.HeatCapacity ~ t_electron_heatcapacity(mp,sim),
                    dTel.ElecPhon ~ t_electron_phononcoupling(mp,sim),
                    dTel.Spatial ~ 0.0,
                    dTel.Source ~laser,
                    dTel.Tel ~ Tel]
        compose(ODESystem(connections,t;name),dTel)
    elseif sim.Systems.NonEqElectrons==true
        egl = length(mp.egrid)
        @named dTel = athem_elec_template(mp.DOS)
        connections =[dTel.Δu ~ electronelectron_internalenergy(mp.DOS,egl)]

        compose(ODESystem(connections,t;name),dTel)
    end
end
"""``
    This is the template for the electonic temperature ODE. It contains a source term where 
    energy injection functions are defined, a spatial term for thermal conductiivity within the
    electronic system. An electron-phonon coupling term for how the electronic temperature interacts
    with the lattice temperature and then the heat capcaity of the electronic system. All variables
    are then later replaced by their relative function for the current simulation. Any functionality
    that is unwanted during a simulation must be set to 0.0 during setup, not ignored.
"""
function ttm_elec_template(;name)
    @variables Tel(t) Source(t) ElecPhon(t) HeatCapacity(t) Spatial(t)

    eqs = [D(Tel) ~ (Source .+ Spatial .+ ElecPhon)./HeatCapacity]

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

function athem_elec_template(DOS;name)
    @variables Tel(t) Δu(t) Δn(t)
    @parameters kB μ 

    eqs = [D(Tel) ~ 1/(c_T(μ,Tel,DOS,kB)*p_μ(μ,Tel,DOS,kB)-p_T(μ,Tel,DOS,kB)*c_μ(μ,Tel,DOS,kB))*(p_μ(μ,Tel,DOS,kB)*Δu-c_μ(μ,Tel,DOS,kB)*Δn)]

    ODESystem(eqs,t;name)
end

function electronelectron_internalenergy(DOS,egl)
    @variables relax_dis(t)[1:egl]
    @parameters egrid[1:egl] μ u0
    return  relaxedelectron_internalenergy(relax_dis,egrid,DOS,μ,u0)
end

function relaxedelectron_internalenergy(relax_dis::AbstractVector,egrid,DOS::Spline1D,μ::Real,u0::Real)
    dis_spl = get_interpolate(egrid,relax_dis)
    return get_internalenergyspl(μ,dis_spl,DOS,u0)
end
@register_symbolic relaxedelectron_internalenergy(relax_dis::AbstractVector,egrid::AbstractVector,DOS::Spline1D,μ::Num,u0::Num)::Real
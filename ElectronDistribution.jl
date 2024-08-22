"""
    This the base factory function that constructs the non-equilibrium electron ODE. All variables and
    functionality for the on-equilibrium electron ODE should be set up within this function call.
"""
function athem_factory(sim::SimulationSettings,DOS::Spline1D,laser::Num,egl::Int;name)
    if sim.Systems.ElectronTemperature == true
        @variables Tel(t)
    else
        @parameters Tel
    end
    @parameters kB egrid[1:egl] μ
    @named dfneq = athem_template(egl)
    feq = FermiDirac(egrid,μ,Tel,kB)
    connections = [dfneq.Source ~ athem_wrapper(dfneq.fneq,feq,DOS,egl).*laser,
                   dfneq.neqelel ~ electronelectron_wrapper(sim,dfneq.fneq,feq,DOS,egl)]

    compose(ODESystem(connections,t;name),dfneq)
end
"""
    This is the template for the non-equilibrium electron ODE. It contains a source term where 
    energy injection functions are defined, as well as terms for electorn-electron and electron-
    phonon interactions. All variables are then later replaced by their relative function for 
    the current simulation. Any functionality that is unwanted during a simulation must be set 
    to 0.0 during setup, not ignored.
"""
function athem_template(egl::Int;name)
    @variables (fneq(t))[1:egl] (Source(t))[1:egl] (neqelel(t))[1:egl] #(neqelph(t))[1:egl]

    eqs = D(fneq) ~ Source .+ neqelel #.+ neqelph

    ODESystem(eqs,t;name)
end

function athem_wrapper(fneq,feq,DOS,egl)
    @parameters (egrid)[1:egl] hv
    return athem_excitation(egrid,feq,fneq,hv,DOS)
end

function athem_excitation(egrid,feq,fneq,hv::Real,DOS::Spline1D)
    ftot = fneq.+feq
    ftotspl=get_interpolate(egrid,ftot)

    Δfe = fgr_electron_generation(egrid,DOS,ftotspl,hv)
    Δfh = fgr_hole_generation(egrid,DOS,ftotspl,hv)

    pc_sf = get_noparticles(Δfe,DOS,egrid) / get_noparticles(Δfh,DOS,egrid)
    Δfshape = (Δfe.*pc_sf).-Δfh
    inten = get_internalenergy_grid(Δfshape,DOS,egrid)

    return Δfshape./inten
end
@register_array_symbolic athem_excitation(egrid::AbstractVector,feq::AbstractVector,fneq::AbstractVector,hv::Num,DOS::Spline1D) begin
    size = (length(egrid),)
    eltype = eltype(fneq)
end


function fgr_hole_generation(egrid,DOS::Spline1D,ftotspl::Spline1D,hv::Real)
    return DOS(egrid.+hv).*ftotspl(egrid).*(1 .-ftotspl(egrid.+hv))
end

function fgr_electron_generation(egrid,DOS::Spline1D,ftotspl::Spline1D,hv::Real)
    return DOS(egrid.-hv).*ftotspl(egrid.-hv).*(1 .-ftotspl(egrid))
end

function electronelectron_wrapper(sim,fneq,feq,DOS,egl)
    if sim.Interactions.ElectronElectron == true
        @variables Tel(t) n(t)
        @parameters egrid[1:egl] μ kB u0 FE τ
        return athem_electronelectronrelaxation(fneq,feq,egrid,μ,kB,u0,n,DOS)./flt_relaxation(τ,FE,μ,egrid,kB,Tel)
    else
        return zeros(egl)
    end
end

function athem_electronelectronrelaxation(fneq,feq,egrid,μ::Real,kB::Real,u0::Real,n::Real,DOS::Spline1D)
    ftot = get_interpolate(egrid,feq.+fneq)
    goal = get_internalenergyspl(μ,ftot,DOS,u0)
    frel = find_relaxed_eedistribution(goal,kB,egrid,n,DOS,u0)
    return (fneq .+ frel .- feq)
end
@register_array_symbolic athem_electronelectronrelaxation(fneq::AbstractVector,feq::AbstractVector,egrid::AbstractVector,μ::Num,kB::Num,u0::Num,n::Num,DOS::Spline1D) begin
    size=(length(fneq),)
    eltype = eltype(fneq)
end

function find_relaxed_eedistribution(goal,kB,egrid,n,DOS,u0)
    f(u) = goal - relaxed_internalenergy(u,n,DOS,kB,egrid,u0)
    rel_Tel = solve(ZeroProblem(f,1000.0),Order1();atol=1e-4,rtol=1e-4)
    rel_μ = find_chemicalpotential(n,rel_Tel,0.0,DOS,kB)
    return FermiDirac(egrid,rel_μ,rel_Tel,kB)
end

function relaxed_internalenergy(T,n,DOS,kB,egrid,u0)
    cp = find_chemicalpotential(n,T,0.0,DOS,kB)
    u = get_internalenergyspl(cp,get_interpolate(egrid,FermiDirac(egrid,cp,T,kB)),DOS,u0)
    return u
end

function flt_relaxation(τ,FE,μ,egrid,kB,Tel)
    return τ*(FE+μ)^2 ./((egrid.-(FE+μ)).^2 .+(pi*kB*Tel)^2)
end

function dFDdE(kB::Float64,Tel::Float64,μ::Float64,E::Float64)::Real
    Numer=-exp((E-μ)/(kB*Tel))
    Denom=kB*Tel*(exp((E-μ)/(kB*Tel))+1)^2
    return Numer/Denom
end

function dFDdT(kB::Float64,Tel::Float64,μ::Float64,E::Float64)::Real
    Numer=(E-μ)*exp((E-μ)/(kB*Tel))
    Denom=kB*Tel^2*(exp((E-μ)/(kB*Tel))+1)^2
    return Numer/Denom
end

function dFDdμ(kB::Float64,Tel::Float64,μ::Float64,E::Float64)::Real
    Numer=exp((E-μ)/(kB*Tel))
    Denom=kB*Tel*(exp((E-μ)/(kB*Tel))+1)^2
    return Numer/Denom
end

FermiDirac(E::Real,μ::Union{Real,ForwardDiff.Dual},Tel::Real,kB::Real)= 1/(exp((E-μ)/(kB*Tel))+1)

FermiDirac(egrid,μ::Real,Tel::Real,kB::Real) = 1 ./(exp.((egrid.-μ)./(kB*Tel)).+1)
@register_array_symbolic FermiDirac(egrid::AbstractVector,μ::Num,Tel::Num,kB::Num) begin
    size = (length(egrid),)
    eltype = eltype(egrid)
end

function particle_change(DOS,egl;name)
    @variables n(t)

    eqs = D(n) ~ electronelectron_particlechange(n,DOS,egl)

    ODESystem(eqs,t;name)
end

function electronelectron_particlechange(n,DOS,egl)
    @variables fneq(t)[1:egl] Tel(t)
    @parameters τ FE egrid[1:egl] μ kB u0 n0
    feq = FermiDirac(egrid,μ,Tel,kB)
    ee_dis = athem_electronelectronrelaxation(fneq,feq,egrid,μ,kB,u0,n,DOS)
    relax_dis = ee_dis./flt_relaxation(τ,FE,μ,egrid,kB,Tel)
    return relaxedelectron_particlechange(relax_dis,egrid,μ,DOS,n0) 
end

function relaxedelectron_particlechange(relax_dis,egrid,μ::Real,DOS::Spline1D,n0::Real)
    relax_spl = get_interpolate(egrid,relax_dis)
    return get_noparticlesspl(μ,relax_spl,DOS,n0)
end
@register_symbolic relaxedelectron_particlechange(relax_dis::AbstractVector,egrid::AbstractVector,μ::Num,DOS::Spline1D,n0::Num)
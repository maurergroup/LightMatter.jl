"""
    This the base factory function that constructs the non-equilibrium electron ODE. All variables and
    functionality for the on-equilibrium electron ODE should be set up within this function call.
"""
function athem_factory(DOS::Spline1D,laser::Num,egl::Int;name)
    @variables (ftot(t))[1:egl] = zeros(egl)
    @named dneqFD = athem_template(egl)
    @named source =  athem_excitation(dneqFD.fneq,DOS,laser,egl)
    connections = [dneqFD.Source ~ source.excite]

    compose(ODESystem(connections,t;name),dneqFD,source)
end
"""
    This is the template for the non-equilibrium electron ODE. It contains a source term where 
    energy injection functions are defined, as well as terms for electorn-electron and electron-
    phonon interactions. All variables are then later replaced by their relative function for 
    the current simulation. Any functionality that is unwanted during a simulation must be set 
    to 0.0 during setup, not ignored.
"""
function athem_template(egl::Int;name)
    @variables begin
        Source(t)[1:egl]
        fneq(t)[1:egl]# (neqelel(t))[1:egl] (neqelph(t))[1:egl]
    end
    eqs = D(fneq) ~ Source# .+ neqelel .+ neqelph

    ODESystem(eqs,t;name)
end

function athem_excitation(fneq,DOS::Spline1D,laser::Num,egl::Int;name)
    @variables Δfneqh(t)[1:egl] Δfneqe(t)[1:egl] pc_sf(t) δ(t) ftot(t)[1:egl] excite(t)[1:egl] Δfshape(t)[1:egl]
    @parameters egrid[1:egl] μ kB Tel hv

    eqs = [ftot .~ fneq .+ FermiDirac(egrid,μ,Tel,kB)
            Δfneqh .~ fgr_hole_generation(egrid,DOS,ftot,hv)
            Δfneqe .~ fgr_electron_generation(egrid,DOS,ftot,hv)
            pc_sf ~ fgr_particleconservation(DOS,Δfneqh,Δfneqe,egrid,μ)
            Δfshape .~ (Δfneqe*pc_sf) .- Δfneqh
            δ ~ laser/fgr_excitation_internalenergy(Δfshape,DOS,egrid,μ)
            excite .~ δ*Δfshape]
    
    ODESystem(eqs,t;name)
    
end

function fgr_hole_generation(egrid::Vector{Float64},DOS::Spline1D,ftot::Vector{Float64},hv::Real)
    println("Test")
    ftotspl=get_interpolate(egrid,ftot)
    return DOS.(egrid.+hv).*ftotspl.(egrid).*(1 .-ftotspl.(egrid.+hv))
end
@register_array_symbolic fgr_hole_generation(egrid::AbstractVector,DOS::Spline1D,ftot::AbstractVector,hv::Num) begin
    size = (length(egrid),)
    eltype = eltype(ftot)
end

function fgr_electron_generation(egrid::Vector{Float64},DOS::Spline1D,ftot::Vector{Float64},hv::Real)
    println("Test")
    ftotspl=get_interpolate(egrid,ftot)
    return DOS.(egrid.-hv).*ftotspl.(egrid.-hv).*(1 .-ftotspl.(egrid))
end
@register_array_symbolic fgr_electron_generation(egrid::AbstractVector,DOS::Spline1D,ftot::AbstractVector,hv::Num) begin
    size = (length(egrid),)
    eltype = eltype(ftot)
end

function fgr_particleconservation(DOS::Spline1D,fneqh::Vector{Float64},fneqe::Vector{Float64},egrid::Vector{Float64},μ::Real)
    println("Test")
    elDis = get_interpolate(egrid,fneqe)
    hDis = get_interpolate(egrid,fneqh)
    f(u,p) = get_noparticles(μ,hDis,DOS) - u*get_noparticles(μ,elDis,DOS)
    return solve(NonlinearProblem(f,1.0),SimpleKlement();abstol=1e-3,reltol=1e-3).u
end
@register_symbolic fgr_particleconservation(DOS::Spline1D,fneqh::AbstractVector,fneqe::AbstractVector,egrid::AbstractVector,μ::Num)::Real

function fgr_excitation_internalenergy(Δfneq::Vector{Real},DOS::Spline1D,egrid::Vector{Float64},μ::Real)
    fneqspl = get_interpolate(egrid,fneq)
    return get_internalenergy(μ,fneqspl,DOS)
end
@register_symbolic fgr_excitation_internalenergy(Δfneq::AbstractVector,DOS::Spline1D,egrid::AbstractVector,μ::Num)

function athem_electronelectron(sim::SimulationSettings,egl::Int)
    if sim.Interactions.ElectronElectron == true
        @variables (neqFD(t))[1:egl]
        @parameters egrid[1:egl] μ Tel kB 
        eqFD = FermiDirac(egrid,μ,Tel,kB)
        FDrel = find_relaxed_eedistriubtion(egl::Int)
        τ = flt_relaxation(egl::Int)
    else 
        return 0.0
    end
end

function find_relaxed_eedistribution()
    return
end

function flt_relaxation(egl::Int)
    @parameters FE τ μ Tel kB egrid[1:egl]
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

FermiDirac(E::Real,μ::Union{Real,ForwardDiff.Dual},Tel::Real,kB::Real)= 1/(exp((E-μ)/(kB*Tel))+1)
FermiDirac(egrid,μ::Num,Tel::Num,kB::Num)= 1 ./(exp.((egrid.-μ)./(kB*Tel)).+1)

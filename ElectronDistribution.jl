"""
    This the base factory function that constructs the non-equilibrium electron ODE. All variables and
    functionality for the on-equilibrium electron ODE should be set up within this function call.
"""
function athem_factory(DOS::Spline1D,laser::Num,egl::Int;name)
    @variables (fneq(t))[1:egl]

    @named dfneq = athem_template(egl)
    
    connections = [dfneq.Source ~ athem_excitation_wrapper(DOS,egl,laser),
                   dfneq.fneq ~ fneq]

    #connections = Symbolics.scalarize.(reduce(vcat,Symbolics.scalarize.(connections)))

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
    @variables (fneq(t))[1:egl] (Source(t))[1:egl]#(neqelel(t))[1:egl] (neqelph(t))[1:egl]

    eqs = D(fneq) ~ Source# .+ neqelel .+ neqelph

    ODESystem(eqs,t;name)
end

function athem_excitation_wrapper(DOS,egl,laser)
    @variables (fneq(t))[1:egl]
    @parameters (egrid)[1:egl] μ Tel kB hv
    return athem_excitation(egrid,μ,Tel,kB,fneq,hv,DOS).*laser
end

function athem_excitation(egrid::AbstractVector,μ::Real,Tel::Real,kB::Real,fneq::AbstractVector,hv::Real,DOS::Spline1D)

    ftot = fneq.+FermiDirac(egrid,μ,Tel,kB)
    ftotspl=get_interpolate(egrid,ftot)

    Δfe = fgr_electron_generation(egrid,DOS,ftotspl,hv)
    Δfh = fgr_hole_generation(egrid,DOS,ftotspl,hv)

    pc_sf = fgr_particleconservation(DOS,Δfh,Δfe,egrid,μ)
    Δfshape = (Δfe.*pc_sf).-Δfh
    inten = fgr_excitation_internalenergy(Δfshape,DOS,egrid,μ)

    return Δfshape./inten
end
@register_array_symbolic athem_excitation(egrid::AbstractVector,μ::Num,Tel::Num,kB::Num,fneq::AbstractVector,hv::Num,DOS::Spline1D) begin
    size = (length(egrid),)
    eltype = eltype(fneq)
end


function fgr_hole_generation(egrid::AbstractVector,DOS::Spline1D,ftotspl::Spline1D,hv::Real)
    return DOS.(egrid.+hv).*ftotspl.(egrid).*(1 .-ftotspl.(egrid.+hv))
end

function fgr_electron_generation(egrid::AbstractVector,DOS::Spline1D,ftotspl::Spline1D,hv::Real)
    return DOS.(egrid.-hv).*ftotspl.(egrid.-hv).*(1 .-ftotspl.(egrid))
end

function fgr_particleconservation(DOS::Spline1D,fneqh::AbstractVector,fneqe::AbstractVector,egrid::AbstractVector,μ::Real)
    elDis = get_interpolate(egrid,fneqe)
    hDis = get_interpolate(egrid,fneqh)
    f(u,p) = get_noparticles(μ,hDis,DOS) - u*get_noparticles(μ,elDis,DOS)
    return solve(NonlinearProblem(f,1.0),SimpleKlement();abstol=1e-3,reltol=1e-3).u
end

function fgr_excitation_internalenergy(Δfneq::AbstractVector,DOS::Spline1D,egrid::AbstractVector,μ::Real)
    fneqspl = get_interpolate(egrid,Δfneq)
    return get_internalenergy(μ,fneqspl,DOS)
end

#= function athem_electronelectron(sim::SimulationSettings,egl::Int)
    if sim.Interactions.ElectronElectron == true
        @variables (neqFD(t))[1:egl]
        @parameters egrid[1:egl] μ Tel kB 
        eqFD = FermiDirac(egrid,μ,Tel,kB)
        FDrel = find_relaxed_eedistriubtion(egl::Int)
        τ = flt_relaxation(egl::Int)
    else return 0.0
    end
end

function find_relaxed_eedistribution()
    return
end

function flt_relaxation(egl::Int)
    @parameters FE τ μ Tel kB egrid[1:egl]
    return τ*(FE+μ)^2 ./((egrid.-(FE+μ)).^2 .+(pi*kB*Tel)^2)
end =#

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

FermiDirac(egrid,μ,Tel,kB)= 1 ./(exp.((egrid.-μ)./(kB*Tel)).+1)
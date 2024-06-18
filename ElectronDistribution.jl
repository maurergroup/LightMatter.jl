"""
    This the base factory function that constructs the non-equilibrium electron ODE. All variables and
    functionality for the on-equilibrium electron ODE should be set up within this function call.
"""
function athem_factory(DOS::Spline1D,laser::Num,egl::Int;name)
    @named dneqFD = athem_template(egl)
    #= connections = [dneqFD.Source ~ athem_excitation(DOS,laser,egl),
                   dneqFD.neqelel ~ athem_electronelectron(sim,egl),
                   dneqFD.neqelph ~ athem_electronphonon(sim,egl)] =#
    connections = [dneqFD.Source ~ athem_excitation(DOS,laser,egl)]
    compose(ODESystem(connections,t;name),dneqFD)
end
"""
    This is the template for the non-equilibrium electron ODE. It contains a source term where 
    energy injection functions are defined, as well as terms for electorn-electron and electron-
    phonon interactions. All variables are then later replaced by their relative function for 
    the current simulation. Any functionality that is unwanted during a simulation must be set 
    to 0.0 during setup, not ignored.
"""
function athem_template(egl::Int;name)
    @variables (neqFD(t))[1:egl] (Source(t))[1:egl] (neqelel(t))[1:egl] (neqelph(t))[1:egl]

    eqs = D(neqFD) ~ Source# .+ neqelel .+ neqelph

    ODESystem(eqs,t;name)
end

function athem_excitation(DOS::Spline1D,laser::Num,egl::Int)
    @variables (Δfneqtot(t))[1:egl] δ(t) (Δfneqe(t))[1:egl] (Δfneqh(t))[1:egl]
    Δfneqh ~ neqhole_generation(DOS,egl)
    Δfneqe ~ athem_electron_shape(DOS,egl)
    δ ~ athem_excitation_size(DOS,laser,egl)
    Δfneqtot ~ δ*(Δfneqe-Δfneqh)
    return Δfneqtot
end
    
function athem_electron_shape(DOS::Spline1D,egl::Int)
    @variables (Δfneqe(t))[1:egl] (Δfneqh(t))[1:egl]
    @parameters egrid[1:egl]

    neqelectron_generation(DOS,egl).*athem_particleconservation(DOS,neqelectron_generation(DOS,egl),Δfneqh,egrid)
end

function athem_excitation_size(DOS::Spline1D,laser::Num,egl::Int)
    @variables (Δfneqe(t))[1:egl] (Δfneqh(t))[1:egl]
    @parameters μ egrid[1:egl]
    Δfneqspl = get_interpolate(egrid,Δfneqe-Δfneqh)
    return laser / get_internalenergy(μ,Δfneqspl,DOS)
end

function neqelectron_generation(DOS::Spline1D,egl::Int)
    @parameters hv μ Tel kB egrid[1:egl]
    DOS.(egrid.-hv).*FermiDirac(egrid.-hv,μ,Tel,kB).*(1 .-FermiDirac(egrid,μ,Tel,kB))
end

function neqhole_generation(DOS::Spline1D,egl::Int)
    @parameters hv μ Tel kB egrid[1:egl]
    DOS.(egrid.+hv).*FermiDirac(egrid,μ,Tel,kB).*(1 .-FermiDirac(egrid.+hv,μ,Tel,kB))
end

function athem_particleconservation(DOS::Spline1D,fneqe::Vector{Float64},fneqh::Vector{Float64},egrid::Vector{Float64})
    elDis = get_interpolate(egrid,fneqe)
    hDis = get_interpolate(egrid,fneqh)
    f(u,p) = get_noparticles(μ,hDis,DOS) - u*get_noparticles(μ,elDis,DOS)
    return solve(NonlinearProblem(f,1.0),SimpleKlement();abstol=1e-3,reltol=1e-3).u
end
@register_array_symbolic athem_particleconservation(DOS::Spline1D,fneqe::Vector{Num},fneqh::Vector{Num},egrid::Vector{Num})::Real begin
    size=1
    eltype=eltype(fneqe)
end

function athem_electronelectron(sim::SimulationSettings,egl::Int)
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
FermiDirac(E::Symbolics.Arr{Num,1},μ::Num,Tel::Num,kB::Num) = 1 ./(exp.((E.-μ)./(kB*Tel)).+1)

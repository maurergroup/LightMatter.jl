include("SimulationVariables.jl")


function athem_factory(DOS::Spline1D,ERange::Vector{<:Real},laser::Num)
    @parameters totalDis(t)
    @named dneqFD = athem_template()
    connections = [dneqFD.Source ~ athem_excitation(DOS,ERange,laser)]
    connected = compose(ODESystem(connections,t;name=:connected),dneqFD)
    return connected
end
function athem_template(;name)
    @variables neqFD(t) Source(t)# neqelel(t) neqelph(t)

    eqs = D(neqFD) ~ Source# .+ neqelel .+ neqelph

    ODESystem(eqs,t;name)
end

function athem_excitation(DOS::Spline1D,ERange::Vector{<:Real},laser::Num)
    @variables Δfneqe Δfneqh Δfneqlas δ
    #= ParentScope(Δfneqe)
    ParentScope(Δfneqe)
    ParentScope(Δfneqe)
    ParentScope(δ) =#

    athem_excitation_shape(DOS,ERange)
    athem_excitation_size(DOS,ERange,laser)
    Δfneqlas = δ*(Δfneqe - Δfneqh)
end
    
function athem_excitation_shape(DOS::Spline1D,ERange::Vector{<:Real})
    @variables Δfneqe Δfneqh Δfneqlas
    #= ParentScope(ParentScope(Δfneqe))
    ParentScope(ParentScope(Δfneqh))
    ParentScope(ParentScope(Δfneqlas)) =#

    Δfneqe = neqelectron_generation(DOS,ERange)
    Δfneqh = neqhole_generation(DOS,ERange)
    scale = athem_particleconservation(DOS,ERange,Δfneqe,Δfneqh)
    Δfneqe = Δfneqe*scale
    Δfneqt = Δfneqe .- Δfneqh
end

function athem_excitation_size(DOS::Spline1D,ERange::Vector{<:Real},laser::Num)
    @variables δ Δfneqlas
    @parameters μ
    #= ParentScope(ParentScope(Δfneqt))
    ParentScope(ParentScope(δ))
    ParentScope(ParentScope(μ)) =#

    Δfneqspl = get_interpolate(ERange,Δfneqlas)
    δ = laser / get_internalenergy(μ,Δfneqspl,DOS)
end

function neqelectron_generation(DOS::Spline1D,ERange::Vector{<:Real})
    @parameters hv μ Tel kB
    return DOS.(ERange.-hv).*FermiDirac.(ERange.-hv,μ,Tel,kB).*(1 .-FermiDirac.(ERange,μ,Tel,kB))
end

function neqhole_generation(DOS::Spline1D,ERange::Vector{<:Real})
    @parameters hv μ Tel kB
    return DOS.(ERange.+hv).*FermiDirac.(ERange,μ,Tel,kB).*(1 .-FermiDirac.(ERange.+hv,μ,Tel,kB))
end

function athem_particleconsveration(DOS::Spline1D,ERange::Vector{<:Real},fneqe::Vector{Real},fneqh::Vector{Real})
    elDis = get_interpolate(ERange,fneqe)
    hDis = get_interpolate(ERange,fneqh)
    f(u,p) = get_noparticles(μ,hDis,DOS) - u*get_noparticles(μ,elDis,DOS)
    return solve(NonlinearProblem(f,1.0),SimpleKlement();abstol=1e-3,reltol=1e-3).u
end
@register_symbolic  athem_particleconsveration(DOS::Spline1D,ERange::Vector{<:Real},fneqe::Vector{Num},fneqh::Vector{Num})

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

FermiDirac(E::Float64,μ::Float64,Tel::Float64,kB::Float64) = 1/(exp((E-μ)/(kB*Tel))+1)

FermiDirac(E::Float64,μ::ForwardDiff.Dual,Tel::Float64,kB::Float64) = 1/(exp((E-μ)/(kB*Tel))+1)
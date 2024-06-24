function athemdistribution_factory(sim::SimulationSettings,laser)
    feq = :(1 ./(exp.((mp.egrid.-μ)./(cons.kB*Tel)).+1))
    ftot = :($feq.+fneq)
    athemexcite = athemexcitation(ftot,laser)
    Elecelec = athem_electronelectroninteraction(sim,ftot)
    Elecphon = athem_electronphononinteraction(sim,ftot)
    return build_athemdistribution(athemexcite,Elecelec,Elecphon)
end

function build_athemdistribution(athemexcite,Elecelec,Elecphon)
    return Expr(:call,:.+,athemexcite,Elecelec,Elecphon)
end

function athemexcitation(ftot::Expr,laser::Expr)
    Δfneqe,Δfneqh = athem_neqelectronandhole(ftot::Expr)
    deltas=(Δfneqe,Δfneqh)
    pc_sf = :(athem_particleconservation(mp.DOS,$Δfneqe,$Δfneqh,mp.egrid,μ))
    Δfneqtot = Expr(:call,:.-,Expr(:call,:*,pc_sf,Δfneqe),Δfneqh)
    δ = :($laser/athem_excitation_internalenergy($Δfneqtot,mp.DOS,mp.egrid,μ))
    return Expr(:call,:*,δ,Δfneqtot)
end

function athem_neqelectronandhole(ftot::Expr)
    ftotspl = :(get_interpolate(mp.egrid,$ftot))
    Δfneqh = :(athem_holegeneration(mp.egrid,mp.DOS,$ftotspl,las.hv))
    Δfneqe = :(athem_electrongeneration(mp.egrid,mp.DOS,$ftotspl,las.hv))
    return Δfneqe,Δfneqh
end

function athem_holegeneration(egrid::Vector{Float64},DOS::Spline1D,ftotspl::Spline1D,hv)
    return DOS.(egrid.+hv).*ftotspl.(egrid).*(1 .-ftotspl.(egrid.+hv))
end

function athem_electrongeneration(egrid::Vector{Float64},DOS::Spline1D,ftotspl::Spline1D,hv)
    return DOS.(egrid.-hv).*ftotspl.(egrid.-hv).*(1 .-ftotspl.(egrid))
end

function athem_particleconservation(DOS::Spline1D,Δfneqe,Δfneqh,egrid::Vector{Float64},μ::Float64)
    elDis = get_interpolate(egrid,Δfneqe)
    hDis = get_interpolate(egrid,Δfneqh)
    f(u,p) = get_noparticles(μ,hDis,DOS) - u*get_noparticles(μ,elDis,DOS)
    return solve(NonlinearProblem(f,1.0),SimpleKlement();abstol=1e-3,reltol=1e-3).u
end

function athem_excitation_internalenergy(Δfneqtot::Vector{Float64},DOS::Spline1D,egrid::Vector{Float64},μ::Real)
    fneqspl = get_interpolate(egrid,Δfneqtot)
    return get_internalenergyspl(μ,fneqspl,DOS)
end

function athem_electronelectroninteraction(sim::SimulationSettings,ftot)
    if sim.Interactions.ElectronElectron == true 
        return athem_electronelectronscattering(ftot)
    else
        return 0.0
    end
end

function athem_electronelectronscattering(ftot)
    ftotspl = :(get_interpolate(mp.egrid,$ftot))
    τee = :(τ*mp.FE^2 ./((mp.egrid.-μ).^2 .+(pi*cons.kB*Tel)^2))
    frel = :(find_relaxeddistribution(mp.egrid,get_internalenergyspl(μ,$ftotspl,mp.DOS),no_part,mp.DOS,cons.kB))
    Δfee = Expr(:call,:-,ftot,frel)
    return Expr(:call,:/,Δfee,τee)
end

function athem_electronphononinteraction(sim::SimulationSettings,ftot)
    if sim.Interactions.ElectronPhonon == true 
        return athem_electronphononscattering(ftot)
    else
        return 0.0
    end
end

function athem_electronphononscattering(ftot)
    τep = :(mp.γ*Tel/mp.g)
    frel = :(:(1 ./(exp.((mp.egrid.-μ)./(cons.kB*Tph)).+1)))
    Δfee = Expr(:call,:-,ftot,frel)
    return Expr(:call,:/,Δfee,τep)
end

function find_relaxeddistribution(egrid::Vector{Float64},goal::Real,no_part::Real,DOS::Spline1D,kB::Real)
    f(u,p) = goal - find_temperatureandμ(u,no_part,DOS,kB)
    Temp = solve(NonlinearProblem(f,1000.0),SimpleKlement();abstol=1e-3,reltol=1e-3).u
    μ = find_chemicalpotential(no_part,Temp,0.0,DOS,kB)
    return FermiDirac.(Temp,μ,kB,egrid)
end

function find_temperatureandμ(Tel::Real,no_part::Real,DOS::Spline1D,kB::Real)
    μ = find_chemicalpotential(no_part,Tel,0.0,DOS,kB)
    return get_internalenergy(μ,Tel,DOS,kB)
end

function neqelectron_electrontransfer()
    feq = :(1 ./(exp.((mp.egrid.-μ)./(cons.kB*Tel)).+1))
    ftot = :($feq.+fneq)
    dis = athem_electronelectronscattering(ftot)
    disspl = :(get_interpolate(mp.egrid,$dis))
    return :(get_internalenergyspl(μ,$disspl,mp.DOS))
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

FermiDirac(Tel,μ,kB,E) = 1/(exp((E-μ)/(kB*Tel))+1)
function athemdistribution_factory(sim::SimulationSettings,laser)
    feq = :(1 ./(exp.((mp.egrid.-μ)./(cons.kB*Tel)).+1))
    ftot = :($feq.+fneq)
    athemexcite = athemexcitation(ftot,laser)
    Elecelec = athem_electronelectroninteraction(sim,feq)
    Elecphon = athem_electronphononinteraction(sim,feq)
    return build_athemdistribution(athemexcite,Elecelec,Elecphon)
end

function build_athemdistribution(athemexcite,Elecelec,Elecphon)
    return Expr(:call,:.+,athemexcite,Elecelec,Elecphon)
end

function athemexcitation(ftot::Expr,laser::Expr)
    Δfneqe,Δfneqh = athem_neqelectronandhole(ftot::Expr)
    pc_sf = :(athem_particleconservation(mp.DOS,$Δfneqe,$Δfneqh,mp.egrid,μ))
    Δfneqtot = Expr(:call,:.-,Expr(:call,:*,pc_sf,Δfneqe),Δfneqh)
    δ = :($laser/athem_excitation_internalenergy($Δfneqtot,mp.DOS,mp.egrid,μ))
    return Expr(:call,:*,δ,Δfneqtot)
end

function athem_neqelectronandhole(ftot::Expr)
    ftotspl = :(get_interpolate(mp.egrid,$ftot))
    Δfneqh = :(athem_holegeneration(mp.egrid,mp.DOS,$ftotspl,las.hv,μ))
    Δfneqe = :(athem_electrongeneration(mp.egrid,mp.DOS,$ftotspl,las.hv,μ))
    return Δfneqe,Δfneqh
end

function athem_holegeneration(egrid::Vector{Float64},DOS::Spline1D,ftotspl::Spline1D,hv,μ::Real)
    return DOS.(egrid.+hv).*ftotspl.(egrid).*(1 .-ftotspl.(egrid.+hv))./DOS(μ)
end

function athem_electrongeneration(egrid::Vector{Float64},DOS::Spline1D,ftotspl::Spline1D,hv,μ::Real)
    return DOS.(egrid.-hv).*ftotspl.(egrid.-hv).*(1 .-ftotspl.(egrid))./DOS(μ)
end

function athem_particleconservation(DOS::Spline1D,Δfneqe,Δfneqh,egrid::Vector{Float64},μ::Float64)
    elDis = get_interpolate(egrid,Δfneqe)
    hDis = get_interpolate(egrid,Δfneqh)
    f(u,p) = get_noparticles(μ,hDis,DOS) - u*get_noparticles(μ,elDis,DOS)
    return solve(NonlinearProblem(f,1.0),Klement();abstol=1e-3,reltol=1e-3).u
end

function athem_excitation_internalenergy(Δfneqtot::Vector{Float64},DOS::Spline1D,egrid::Vector{Float64},μ::Real)
    fneqspl = get_interpolate(egrid,Δfneqtot)
    return get_internalenergyspl(μ,fneqspl,DOS)
end

function athem_electronelectroninteraction(sim::SimulationSettings,feq)
    if sim.Interactions.ElectronElectron == true 
        return athem_electronelectronscattering(feq)
    else
        return 0.0
    end
end

function athem_electronelectronscattering(feq)
    ftotspl = :(get_interpolate(mp.egrid,$feq.+fneq))
    τee = :(mp.τ*μ^2 ./((mp.egrid.-μ).^2 .+(pi*cons.kB*Tel)^2))
    frel = :(find_relaxeddistribution(mp.egrid,get_internalenergyspl(μ,$ftotspl,mp.DOS),noe,mp.DOS,cons.kB,mp.FE))
    Δfee = Expr(:call,:.-,frel,:(fneq.+$feq))
    return Expr(:call,:./,Δfee,τee)
end

function athem_electronelectronscattering()
    τee = :(mp.τ*μ^2 ./((mp.egrid.-μ).^2 .+(pi*cons.kB*Tel)^2))
    return Expr(:call,:./,:(fneq*-1),τee)
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
    frel = :(1 ./(exp.((mp.egrid.-μ)./(cons.kB*Tph)).+1))
    Δfee = Expr(:call,:-,ftot,frel)
    return Expr(:call,:/,Δfee,τep)
end

function find_relaxeddistribution(egrid::Vector{<:Real},goal::Real,no_part::Real,DOS::Spline1D,kB::Real,FE::Real)
    f(u) = goal - find_temperatureandμ(u,no_part,DOS,kB,FE)
    Temp = solve(ZeroProblem(f,1000.0),Order1();abstol=1e-6,reltol=1e-6)
    μ = find_chemicalpotential(no_part,Temp,FE,DOS,kB)
    return FermiDirac.(Temp,μ,kB,egrid)
end

function find_temperatureandμ(Tel::Real,no_part::Real,DOS::Spline1D,kB::Real,FE::Real)
    μ = find_chemicalpotential(no_part,Tel,FE,DOS,kB)
    return get_internalenergy(μ,Tel,DOS,kB)
end

function neqelectron_electrontransfer()
    feq = :(1 ./(exp.((mp.egrid.-μ)./(cons.kB*Tel)).+1))
    dis = athem_electronelectronscattering(feq)
    disspl = :(get_interpolate(mp.egrid,$dis))
    inten = :(get_internalenergyspl(μ,$disspl,mp.DOS))
    return inten
end

function neqelectron_electronparticlechange()
    feq = :(1 ./(exp.((mp.egrid.-μ)./(cons.kB*Tel)).+1))
    dis = athem_electronelectronscattering(feq)
    disspl = :(get_interpolate(mp.egrid,$dis))
    return :(get_noparticles(μ,$disspl,mp.DOS))
end

function neqelectron_phonontransfer()
    feq = :(1 ./(exp.((mp.egrid.-μ)./(cons.kB*Tel)).+1))
    ftot = :($feq.+fneq)
    dis = athem_electronphononscattering(ftot)
    disspl = :(get_interpolate(mp.egrid,$dis))
    return :(get_internalenergyspl(μ,$disspl,mp.DOS))
end

function dFDdE(kB::Float64,Tel::Real,μ::Float64,E::Float64)::Real
    Numer=-exp((E-μ)/(kB*Tel))
    Denom=kB*Tel*(exp((E-μ)/(kB*Tel))+1)^2
    return Numer/Denom
end

function dFDdT(kB::Float64,Tel::Real,μ::Float64,E::Float64)::Real
    Numer=(E-μ)*exp((E-μ)/(kB*Tel))
    Denom=kB*Tel^2*(exp((E-μ)/(kB*Tel))+1)^2
    return Numer/Denom
end

FermiDirac(Tel,μ,kB,E) = 1 ./(exp.((E.-μ)./(kB*Tel)).+1)


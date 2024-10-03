function athemdistribution_factory(sim::SimulationSettings,laser,dim)
    feq = :(1 ./(exp.((mp.egrid.-μ)./(cons.kB*Tel)).+1))
    ftot = :($feq.+fneq)
    Elecelec = athem_electronelectroninteraction(sim)
    Elecphon = athem_electronphononinteraction(sim)
    athemexcite=athemexcitation(:(fneq.+$feq),laser)
    return build_athemdistribution(athemexcite,Elecelec,Elecphon)
end

function build_athemdistribution(athemexcite,Elecelec,Elecphon)
    return Expr(:call,:.+,athemexcite,Elecelec,Elecphon)
end

function athemexcitation(ftot::Expr,laser::Expr)
    ftotspl = :(get_interpolate(mp.egrid,$ftot))
    Δfneqh = :(athem_holegeneration(mp.egrid,mp.DOS,$ftotspl,las.hv))
    Δfneqe = :(athem_electrongeneration(mp.egrid,mp.DOS,$ftotspl,las.hv))
    pc_sf = :(athem_particleconservation(mp.DOS,$Δfneqe,$Δfneqh,mp.egrid,μ,mp.n0,mp.FE))
    Δfneqtot = Expr(:call,:.-,Expr(:call,:*,pc_sf,Δfneqe),Δfneqh)
    δ = :($laser/athem_excitation_internalenergy($Δfneqtot,mp.DOS,mp.egrid,μ,mp.u0,mp.FE))
    return Expr(:call,:*,δ,Δfneqtot)
end

function athem_holegeneration(egrid::Vector{Float64},DOS::spl,ftotspl,hv)
    return DOS(egrid.+hv).*ftotspl(egrid).*(1 .-ftotspl(egrid.+hv))
end

function athem_electrongeneration(egrid::Vector{Float64},DOS::spl,ftotspl,hv)
    return DOS(egrid.-hv).*ftotspl(egrid.-hv).*(1 .-ftotspl(egrid))
end

function athem_particleconservation(DOS::spl,Δfneqe,Δfneqh,egrid::Vector{Float64},μ::Float64,n0,FE)
    elDis = get_interpolate(egrid,Δfneqe)
    hDis = get_interpolate(egrid,Δfneqh)
    f(u) = get_noparticlesspl(μ,hDis,DOS,n0,FE) - u*get_noparticlesspl(μ,elDis,DOS,n0,FE)
    return solve(ZeroProblem(f,1.0),Order1();atol=1e-2,rtol=1e-2)
end

function athem_excitation_internalenergy(Δfneqtot::Vector{Float64},DOS::spl,egrid::Vector{Float64},μ::Real,u0,FE)
    fneqspl = get_interpolate(egrid,Δfneqtot)
    return get_internalenergyspl(μ,fneqspl,DOS,u0,FE)
end

function athem_electronelectroninteraction(sim::SimulationSettings)
    if sim.Interactions.ElectronElectron == true 
        return :(-1*relax_dis)
    else
        return 0.0
    end
end

function athem_electronelectronscattering()
    feq = :(FermiDirac(Tel,μ,cons.kB,mp.egrid))
    ftotspl = :(get_interpolate(mp.egrid,$feq.+fneq))
    τee = :(mp.τ*(μ.+mp.FE)^2 ./((mp.egrid.-μ).^2 .+(pi*cons.kB*Tel)^2))
    goal = :(get_internalenergyspl(μ,$ftotspl,mp.DOS,mp.u0,mp.FE))
    frel = :(find_relaxeddistribution(mp.egrid,$goal,n,mp.DOS,cons.kB,mp.u0,μ,mp.FE))
    Δfee = Expr(:call,:.-,:(fneq.+$frel),feq)
    return Expr(:call,:./,Δfee,τee)
end

function find_relaxeddistribution(egrid::Vector{<:Real},goal::Real,n::Real,DOS::spl,kB::Real,u0,μ,FE)
    f(u) = goal - find_temperatureandμ(u,n,DOS,kB,u0,μ,egrid,FE)
    Temp = solve(ZeroProblem(f,1000.0),Order1();abstol=1e-10,reltol=1e-10)
    μ = find_chemicalpotential(n,Temp,DOS,kB,FE)
    return FermiDirac(Temp,μ,kB,egrid)
end

function find_temperatureandμ(Tel::Real,n::Real,DOS::spl,kB::Real,u0::Real,μ,egrid,FE)
    spl = get_interpolate(egrid,FermiDirac(Tel,μ,kB,egrid))
    μ = find_chemicalpotential(n,Tel,DOS,kB,FE)
    return get_internalenergyspl(μ,spl,DOS,u0,FE)
end

function athem_electronphononinteraction(sim::SimulationSettings)
    if sim.Interactions.ElectronPhonon == true 
        return :(-fneq./mp.τep)
    else
        return 0.0
    end
end

function athem_electronparticlechange()
    spl = :(get_interpolate(mp.egrid,relax_dis))
    return :(get_noparticlesspl(μ,$spl,mp.DOS,mp.n0,mp.FE))
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

function dFDdμ(kB::Float64,Tel::Float64,μ::Float64,E::Float64)::Real
    Numer=exp((E-μ)/(kB*Tel))
    Denom=kB*Tel*(exp((E-μ)/(kB*Tel))+1)^2
    return Numer/Denom
end

FermiDirac(Tel,μ,kB,E) = 1 ./(exp.((E.-μ)./(kB*Tel)).+1)


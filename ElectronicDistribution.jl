function athemdistribution_factory(sim::SimulationSettings,laser::Expr)
    feq = :(1 ./(exp.((mp.egrid.-μ)./(cons.kB*Tel)).+1))
    ftot = :($feq.+fneq)
    Elecelec = athem_electronelectroninteraction(sim)
    Elecphon = athem_electronphononinteraction(sim)
    athemexcite=:($laser*(athemexcitation($ftot,mp.egrid,mp.DOS,las.hv,mp.n0,mp.FE,mp.u0)))
    return build_athemdistribution(athemexcite,Elecelec,Elecphon)
end

function build_athemdistribution(athemexcite,Elecelec,Elecphon)
    return Expr(:call,:.+,athemexcite,Elecelec,Elecphon)
end

function athemexcitation(ftot,egrid,DOS,hv,n0,FE,u0)
    ftotspl = get_interpolate(egrid,ftot)
    Δfneqh = athem_holegeneration(egrid,DOS,ftotspl,hv)
    Δfneqe = athem_electrongeneration(egrid,DOS,ftotspl,hv)
    pc_sf = get_noparticlesspl(get_interpolate(egrid,Δfneqe),DOS,n0,FE) / get_noparticlesspl(get_interpolate(egrid,Δfneqh),DOS,n0,FE)
    Δfneqtot = (pc_sf*Δfneqe).-Δfneqh
    return Δfneqtot./get_internalenergyspl(get_interpolate(egrid,Δfneqtot),DOS,u0,FE)
end

function athem_holegeneration(egrid::Vector{Float64},DOS::spl,ftotspl::spl,hv::Float64)
    return DOS(egrid.+hv).*ftotspl(egrid).*(1 .-ftotspl(egrid.+hv))
end

function athem_electrongeneration(egrid::Vector{Float64},DOS::spl,ftotspl::spl,hv::Float64)
    return DOS(egrid.-hv).*ftotspl(egrid.-hv).*(1 .-ftotspl(egrid))
end

function athem_particleconservation(DOS::spl,Δfneqe::Vector{Float64},Δfneqh::Vector{Float64},egrid::Vector{Float64},n0::Float64,FE::Float64)
    elDis = get_interpolate(egrid,Δfneqe)
    hDis = get_interpolate(egrid,Δfneqh)
    f(u) = get_noparticlesspl(hDis,DOS,n0,FE) - u*get_noparticlesspl(elDis,DOS,n0,FE)
    return solve(ZeroProblem(f,1.0),Order1();atol=1e-2,rtol=1e-2)
end

function athem_electronelectroninteraction(sim::SimulationSettings)
    if sim.Interactions.ElectronElectron == true 
        return :(-1*relax_dis)
    else
        return 0.0
    end
end

function athem_electronphononinteraction(sim::SimulationSettings)
    if sim.Interactions.ElectronElectron == true 
        return :(-1*fneq./mp.τep)
    else
        return 0.0
    end
end

function athem_electronelectronscattering()
    feq = :(FermiDirac(Tel,μ,cons.kB,mp.egrid))
    ftotspl = :(get_interpolate(mp.egrid,$feq.+fneq))
    τee = :(mp.τ*(μ.+mp.FE)^2 ./((mp.egrid.-μ).^2 .+(pi*cons.kB*Tel)^2))
    goal = :(get_internalenergyspl($ftotspl,mp.DOS,mp.u0,mp.FE))
    frel = :(find_relaxeddistribution(mp.egrid,$goal,n,mp.DOS,cons.kB,mp.u0,mp.FE,mp.n0))
    Δfee = Expr(:call,:.-,:(fneq.+$frel),feq)
    return Expr(:call,:./,Δfee,τee)
end

function find_relaxeddistribution(egrid::Vector{Float64},goal::Float64,n::Float64,DOS::spl,kB::Float64,u0::Float64,FE::Float64,n0::Float64)
    f(u) = goal - find_temperatureandμ(u,n,DOS,kB,u0,egrid,FE,n0)
    Temp = solve(ZeroProblem(f,1000.0);abstol=1e-7,reltol=1e-7)
    μ = find_chemicalpotential(n,Temp,DOS,kB,FE,n0)
    return FermiDirac(Temp,μ,kB,egrid)
end

function find_temperatureandμ(Tel::Real,n::Real,DOS::spl,kB::Real,u0::Real,egrid::Vector{Float64},FE::Float64,n0::Float64)
    μ = find_chemicalpotential(n,Tel,DOS,kB,FE,n0)
    spl = get_interpolate(egrid,FermiDirac(Tel,μ,kB,egrid))
    return get_internalenergyspl(spl,DOS,u0,FE)
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
    return :(get_noparticlesspl($spl,mp.DOS,mp.n0,mp.FE))
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

@inline FermiDirac(Tel::Float64,μ::Float64,kB::Float64,E::Union{Vector{Float64},Float64}) = 1 ./(exp.((E.-μ)./(kB*Tel)).+1)


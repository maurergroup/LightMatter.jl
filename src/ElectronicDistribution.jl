function athemdistribution_factory(sim::SimulationSettings,laser::Expr)
    feq = :(1 ./(exp.((mp.egrid.-μ)./(cons.kB*Tel)).+1))
    ftot = :($feq.+fneq)
    Elecelec = athem_electronelectroninteraction(sim)
    Elecphon = athem_electronphononinteraction(sim)
    athemexcite=:($laser*athemexcitation($ftot,mp.egrid,DOS,las.hv))
    return build_athemdistribution(athemexcite,Elecelec,Elecphon)
end

function build_athemdistribution(athemexcite,Elecelec,Elecphon)
    return Expr(:call,:.+,athemexcite,Elecelec,Elecphon)
end

function athemexcitation(ftot,egrid,DOS,hv)
    ftotspl = get_interpolate(egrid,ftot)
    Δfneqh = athem_holegeneration(egrid,DOS,ftotspl,hv)
    Δfneqe = athem_electrongeneration(egrid,DOS,ftotspl,hv)
    pc_sf = get_noparticles(Δfneqe,DOS,egrid) / get_noparticles(Δfneqh,DOS,egrid)
    Δfneqtot = (pc_sf*Δfneqe).-Δfneqh
    return Δfneqtot./get_internalenergy(Δfneqtot,DOS,egrid)
end

function athem_holegeneration(egrid::Vector{Float64},DOS::spl,ftotspl::spl,hv::Float64)
    return DOS(egrid.+hv).*ftotspl(egrid).*(1 .-ftotspl(egrid.+hv))
end

function athem_electrongeneration(egrid::Vector{Float64},DOS::spl,ftotspl::spl,hv::Float64)
    return DOS(egrid.-hv).*ftotspl(egrid.-hv).*(1 .-ftotspl(egrid))
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
    ftot = :($feq.+fneq)
    τee = :(mp.τ*(μ.+mp.FE)^2 ./((mp.egrid.-μ).^2 .+(pi*cons.kB*Tel)^2))
    goal = :(extended_Bode($ftot.*DOS(mp.egrid).*mp.egrid,mp.egrid))
    frel = :(find_relaxeddistribution(mp.egrid,$goal,n,DOS,cons.kB))
    Δfee = Expr(:call,:.-,:(fneq.+$frel),feq)
    return Expr(:call,:./,Δfee,τee)
end

function find_relaxeddistribution(egrid::Vector{Float64},goal::Float64,n::Float64,DOS::spl,kB::Float64)
    f(u) = goal - find_temperatureandμ(u,n,DOS,kB,egrid)
    Temp = solve(ZeroProblem(f,1000.0);abstol=1e-10,reltol=1e-10)
    μ = find_chemicalpotential(n,Temp,DOS,kB,egrid)
    return FermiDirac(Temp,μ,kB,egrid)
end

function find_temperatureandμ(Tel::Real,n::Real,DOS::spl,kB::Real,egrid::Vector{Float64})
    μ = find_chemicalpotential(n,Tel,DOS,kB,egrid)
    return get_internalenergy(FermiDirac(Tel,μ,kB,egrid),DOS,egrid)
end

function athem_electronphononinteraction(sim::SimulationSettings)
    if sim.Interactions.ElectronPhonon == true 
        return :(-fneq./mp.τep)
    else
        return 0.0
    end
end

function athem_electronparticlechange()
    return :([get_noparticles(relax_dis,DOS,mp.egrid)])
end

@inline function dFDdE(kB::Float64,Tel::Real,μ::Float64,E::Union{Vector{Float64},Float64})
    return -exp.((E.-μ)./(kB*Tel))./(kB*Tel*(exp.((E.-μ)./(kB*Tel)).+1).^2)
end

@inline function dFDdT(kB::Float64,Tel::Real,μ::Float64,E::Union{Vector{Float64},Float64})
    return (E.-μ).*exp.((E.-μ)./(kB*Tel))./(kB*Tel^2*(exp.((E.-μ)./(kB*Tel)).+1).^2)
end

@inline function dFDdμ(kB::Float64,Tel::Float64,μ::Float64,E::Union{Vector{Float64},Float64})
    return exp.((E.-μ)./(kB*Tel))./(kB*Tel*(exp.((E.-μ)./(kB*Tel)).+1).^2)
end

@inline FermiDirac(Tel::Float64,μ::Float64,kB::Float64,E::Union{Vector{Float64},Float64}) = 1 ./(exp.((E.-μ)./(kB*Tel)).+1)



function AthEM_Distribution(;name)
    @variables neqFD(t) neqSource(t) neqerelax(t) neqprelax(t)

    eqs = D(neqFD) ~ neqSource .+ neqerelax .+ neqprelax

    ODESystem(eqs,t;name)
end

function dneqdt_laser()
    elecchange=neqelec.(misc.ERange,hv,Tel,μ,gc)
    holechange=neqhole.(misc.ERange,hv,Tel,μ,gc)
    part_scale= neq_particleconservation(elecchange,holechange,mp,misc)
    totalchange=elecchange*part_scale - holechange
    int_en = stepsize(get_interpolate(misc.ERange,totalchange),laser,misc,mp)
    return totalchange*int_en
end

function neq_particleconservation(elecchange,holechange,mp,misc)
    noelec=get_noparticles_int(get_interpolate(misc.ERange,elecchange),mp,misc)
    nohole=get_noparticles_int(get_interpolate(misc.ERange,holechange),mp,misc)
    scale=nohole/noelec
    return scale
end

function neqelec(E,hv,Tel,μ,gc)
    return mp.DOS(E-hv)*FermiDirac(E-hv,μ,Tel,gc)*(1-FermiDirac(E,μ,Tel,gc))/mp.DOS(μ)
end

function neqhole(E,hv,Tel,μ,gc)
    return mp.DOS(E+hv)*FermiDirac(E,μ,Tel,gc)*(1-FermiDirac(E+hv,μ,Tel,gc))/mp.DOS(μ)
end

function stepsize(neqstep,laser,misc,mp)
    int(u,p)=neqstep(u)*mp.DOS(u)*u
    prob=IntegralProblem(int,misc.int_lb,misc.int_ub)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-7,abstol=1e-7)
    return laser/sol.u
end

function neq_elec_relax(neqFD,misc,μ,Tel,gc,mp)
    eqFD = FermiDirac.(misc.ERange,μ,Tel,gc)
    int_en = get_internal_energy_int(get_interpolate(misc.ERange,eqFD.+neqFD),mp,misc)
    relFD= find_relaxed_distribution(int_en,no_part,μ,mp,misc,gc)
    lifetime = FLT_lifetime(E,μ,Tel,gc,τ)
    return -(neqFD+eqFD-relFD)/lifetime
end

function neq_phon_relax(neqFD,misc,μ,Tel,gc,mp,Tph)
    eqFD = FermiDirac.(misc.ERange,μ,Tel,gc)
    relFD= FermiDirac.(misc.ERange,μ,Tph,gc)
    lifetime = neqelec_phon_lifetime(mp,g,Tel,Tph)
    return -(neqFD+eqFD-relFD)/lifetime
end

function find_relaxed_distribution(int_en,no_part,μ,mp,misc,gc)
    f(x) = int_en - relaxed_internal_energy(no_part,x,μ,mp,misc,gc)
    prob=ZeroProblem(f,1000)
    rel_temp = solve(prob,Order16();atol=1e-5,reltol=1e-5)
    rel_μ = find_chemicalpotential(no_part,rel_temp,μ,mp,misc,gc)
    return FermiDirac.(misc.ERange,rel_μ,rel_temp,gc)
end

function relaxed_internal_energy(no_part,Temp,μ,mp,misc,gc)
    y=find_chemicalpotential(no_part,Temp,μ,mp,misc,gc)
    return get_internal_energy(y,Temp,mp,misc)
end

function FLT_lifetime(E,μ,Tel,gc,τ)
    return τ*μ^2 /((E-μ)^2 +(pi*gc.kB*Tel)^2)
end

function FLT_lifetime_prefactor(mp)
    return 128/(sqrt(3)*pi^2*mp.ω)
end

function neqelec_phon_lifetime(mp,g,Tel,Tph)
    return (mp.γ*(Tel+Tph)/(2*g))
end

function neq_elec_int_en(neqFD,misc,μ,Tel,gc,mp)
    neq = neq_elec_relax(neqFD,misc,μ,Tel,gc,mp)
    neqint=get_interpolate(misc.ERange,neq)
    return get_internal_energy_int(neqint,mp,misc)
end

function neq_phonon_int_en(neqFD,misc,μ,Tel,gc,mp,Tph)
    neq = neq_phon_relax(neqFD,misc,μ,Tel,gc,mp,Tph)
    neqint=get_interpolate(misc.ERange,neq)
    return get_internal_energy_int(neqint,mp,misc)
end
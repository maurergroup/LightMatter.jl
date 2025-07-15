
function dynamical_system(du, u, p, t, type::Val{:AthEM})
    athem_conductivity!(p, u)
    fneq, Tel, Tph, n = extract_variables(u, p, Val(p.sim.method))
    athem_loop(du, u, p.sim, t, sim.structure.Elemental_System ≥ 2)
end

function athem_conductivity!(p, u)
    electronictemperature_conductivity!(p, u, Val(sim.electronictemperature.Conductivity))
    phononictemperature_conductivity!(p, u, Val(sim.phononictemperature.Conductivity))
    electronicdistribution_conductivity!(p, u, Val(sim.athermalelectrons.Conductivity), Val(sim.method))
end

function athem_loop(du, u, sim, t, complex::Val{false})
    Threads.@threads for i in eachindex(sim.structure.dimension.grid)
        sim = sim
        DOS = get_DOS(sim, )
        μ = find_chemicalpotential(u.noe, u.Tel, DOS, egrid, sim.ChemicalPotential)
end

function athem_loop(du, u, p, t, complex::Val{true})
    μ = find_chemicalpotential(u.noe, u.Tel, DOS, egrid, sim.ChemicalPotential)
end

function extract_variables(u, p, type::Val{:AthEM})
    fneq = u.fneq
    Tel = 
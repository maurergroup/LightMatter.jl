function boltzmann_E_phononelectron()
    kgrid = sim.stucture.bandstructure[2](sim.structure.egrid)
    κ = boltzmann_screening(f, kgrid, sim)
    dgdt = zeros(length(g))
    for i in eachindex(g)
        Eq = sim.structure.egrid[i]
        q = Eq / (Constants.ħ * sim.phononicdistribution.cs)
        prefac = pi^3 / (3*Constants.ħ*q)
        M = electron_phonon_matrix(sim, Eq, q, κ)
        int(u,p) = boltzmann_pescatter_int(u, q, Eq, sim, f, g, DOS)
        prob = IntegralProblem(int, 0.0, sim.phononicdistribution.ED)
        sol = solve(prob, HCubatureJL(initdiv=10), abstol=1e-3, reltol=1e-3)
        dgdt[i] = sol.u * M * prefac
    end
end

function boltzmann_pescatter_int(E, q, Eq, sim, f, g, DOS)
    tmp = 0.0
    k = sim.structure.bandstructure[2](E)
    for l in -3:3
        J = average_Bessel(l, γ*q)
        D = DOS(E) / k
        for i in [-1, 1]
            E1 = E + (i*Eq) - l*sim.laser.hv
            k1 = sim.structure.bandstructure[2](E1)
            D1 = DOS(E1) / k1
            F = pauli_epscatter_blocking(f, g, i, E, E1)
            Ξ = boltzmann_step(q, k, k1)
            tmp += D * D1 * F * J * Ξ
        end
    end
    return tmp
end
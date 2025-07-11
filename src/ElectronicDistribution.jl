"""
    FermiDirac(Tel::Float64, μ::Float64, E::Union{Vector{Float64},Float64}) 
    
    Returns the thermal occupation of electrons at a given temperature, chemical potential and energy.
    If a vector of energies is given then it will return the distribution across that range.

    # Arguments
    - 'Tel': Thermal electronic temperature
    - 'μ': Chemical Potential
    - 'E': Energy value or range

    # Returns
    - Value or vector of the Fermi-Dirac distribution
"""
@inline FermiDirac(Tel, μ, E) = 1 ./ (exp.((E.-μ) ./ (Constants.kB*Tel)).+1)
"""
    dFDdE(Tel::Float64, μ::Float64, E::Union{Vector{Float64},Float64})
    
    Returns the change in the Fermi distribution with respect to energy at the given energy value or range.

    # Arguments
    - 'Tel': Thermal electronic temperature
    - 'μ': Chemical Potential
    - 'E': Energy value or range

    # Returns
    - Value or vector of the partial derivative of the Fermi distribution with respect to energy
"""
@inline function dFDdE(Tel, μ, E)
    return -exp.((E.-μ)./(Constants.kB*Tel)) ./ (Constants.kB*Tel * (exp.((E.-μ)./(Constants.kB*Tel)).+1).^2)
end
"""
    dFDdT(Tel::Float64, μ::Float64, E::Union{Vector{Float64},Float64})
    
    Returns the change in the Fermi distribution with respect to temperature at the given energy value or range.

    # Arguments
    - 'Tel': Thermal electronic temperature
    - 'μ': Chemical Potential
    - 'E': Energy value or range

    # Returns
    - Value or vector of the partial derivative of the Fermi distribution with respect to temperature
"""
@inline function dFDdT(Tel, μ, E)
    return (E.-μ) .* exp.((E.-μ)./(Constants.kB*Tel)) ./ (Constants.kB*Tel^2 * (exp.((E.-μ)./(Constants.kB*Tel)).+1).^2)
end
"""
    dFDdμ(Tel::Float64, μ::Float64, E::Union{Vector{Float64},Float64})
    
    Returns the change in the Fermi distribution with respect to chemical potential at the given energy value or range.

    # Arguments
    - 'Tel': Thermal electronic temperature
    - 'μ': Chemical Potential
    - 'E': Energy value or range

    # Returns
    - Value or vector of the partial derivative of the Fermi distribution with respect to chemical potential
"""
@inline function dFDdμ(Tel, μ, E)
    return exp.((E.-μ)./(Constants.kB*Tel)) ./ (Constants.kB*Tel * (exp.((E.-μ)./(Constants.kB*Tel)).+1).^2)
end

function magnetotransport_equations(sim)
    B = :($(sim.structure.fields.laser.magnetic) + $(sim.structure.fields.external.magnetic))
    return :(Lightmatter.magnetotransport_1d!(Δf_mt, fneq.+Lightmatter.FermiDirac(Tel, μ, sim.structure.egrid), sim, $B, DOS, n, band, sim.structure.egrid, g_k, tmp))
end

function df_dk!(dfdk::Vector{Float64}, f::Vector{Float64}, bandstructure::Vector{<:AkimaInterpolation}, egrid::Vector{Float64})
    k = bandstructure[2](egrid)
    fspl = DataInterpolations.AkimaInterpolation(f, k, extrapolation = ExtrapolationType.Constant)

    @inbounds for i in eachindex(k)
        dfdk[i] = DataInterpolations.derivative(fspl, k[i])
    end

    return dfdk
end

function magnetotransport_1d!(Δf_mt::Vector{Float64}, f::Vector{Float64}, sim::Simulation, B::Float64, DOS::spl, n::Float64, bandstructure::Vector{<:AkimaInterpolation}, egrid::Vector{Float64}, g_k::Vector{Float64}, dfdk::Vector{Float64})
    h_2_e = get_h2e(sim)

    goal = Lightmatter.get_internalenergy(f, DOS, sim.structure.egrid)
    find_relaxeddistribution(g_k, sim.structure.egrid, goal, n, DOS)
    @inbounds @simd for i in eachindex(f)
        g_k[i] = f[i] - g_k[i]  # g_k now holds f - f₀
    end

    df_dk!(dfdk, g_k, bandstructure, egrid)  # Compute derivative into dfdk (no aliasing)

    v_g = sim.athermalelectrons.v_g
    factor = Constants.q / (Constants.ħ * Constants.c) * B

    @inbounds @simd for i in 1:h_2_e
        Δf_mt[i] = -factor * v_g[i] * dfdk[i] * -1
    end
    @inbounds @simd for i in (h_2_e + 1):length(Δf_mt)
        Δf_mt[i] = factor * v_g[i] * dfdk[i] * -1
    end
end

function get_h2e(sim::Simulation)
    egrid = sim.structure.egrid
    h_2_e = 1
    best = abs(egrid[1])
    @inbounds for i in 2:length(egrid)
        v = abs(egrid[i])
        if v < best
            best = v
            h_2_e = i
        end
    end
    return h_2_e
end

function pauli_blocking(f,E,E_prime)
    return (f(E_prime) * (1-f(E))) - (f(E) * (1-f(E_prime)))
end

function delta_peak(val)
    return isapprox(val, 0.0; rtol=1e-6, atol=1e-6)
end

function boltzmann_E_excitation(f, sim, E_mag)
    kgrid = sim.stucture.bandstructure[2](sim.structure.egrid)
    γ = Bessel_gamma(E_mag, sim)
    κ = Boltzmann_screening(f, kgrid, sim)
    fspl = get_interpolant(sim.structure.egrid, f)
    dfdt = zeros(length(f))
    for i in eachindex(f)
        k = kgrid[i]
        prefac = pi/Constants.ħ / k
        for l in 1:3
            E1 = sim.structure.egrid[i] + l*sim.laser.hv
            k1 = sim.stucture.bandstructure[2](E1)
            dfrac = sim.structure.DOS(E1) / k1
            F = pauli_blocking(fspl, sim.structure.egrid[i], E1)
            mom_change = Boltzmann_E_momentumintegral(l, k, k1, κ, γ, sim)
            dfdt[i] += dfrac * F * mom_change
        end
        dfdt[i] *= prefac
    end
    return dfdt 
end

function Boltzmann_E_momentumintegral(l, k, k1, κ, γ, sim)
    int(u,p) = u* excitation_matrix(sim, u, κ)* average_Bessel(l, γ*u) * Boltzmann_step(u, k, k1)
    prob = IntegralProblem(int, 0.0, 2*k)
    sol = solve(prob, HCubatureJL(initdiv=100), abstol=1e-4, reltol=1e-4)
    return sol.u
end

function electron_electron_matrix(sim, Δk, κ)
    frac1 = Constants.q^2 / Constants.ϵ0 / sim.electronicdistribution.Ω
    frac2 = 1/ (Δk^2 + κ^2)
    return (frac1*frac2)^2
end

function average_Bessel(order, value)
    int(u,p) = besselj(order, value*u)^2
    prob = IntegralProblem(int, -1, 1)
    sol = solve(prob, HCubatureJL(initdiv=50), abstol=1e-5, reltol=1e-5)
    return 1/2 * sol.u
end

function Bessel_gamma(E_mag, sim)
    return Constants.q * E_mag / (sim.electronicdistribution.me * photon_energytofrequency(sim.laser.hv))
end

function Boltzmann_step(q, k, k1)
    val = (q^2 + k^2 - k1^2) / (2*k*q)
    if -1 ≤ val ≤ 1
        return 1.0
    else 
        return 0.0
    end
end

function Boltzmann_screening(f, kgrid, sim)
    prefac = Constants.q^2 * sim.electronicdistribution.me /(pi^2*Constants.ħ^2*Constants.ϵ0)
    return prefac * Bode_rule(f, kgrid)
end
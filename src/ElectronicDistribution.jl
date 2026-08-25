
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
function FermiDirac(Tel, μ, E)
    return 1 ./ (exp.((E.-μ) ./ (Constants.kB*Tel)).+1)
end

function FermiDirac!(dis, Tel, μ, E)
    @. dis = 1 / (exp((E-μ) / (Constants.kB*Tel))+1)
end

function FermiDirac_mul!(dis, Tel, μ, E)
    @. dis = dis * 1 / (exp((E-μ) / (Constants.kB*Tel))+1)
end
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
function dFDdE(Tel, μ, E)
    return -1/(Constants.kB*Tel) .* FermiDirac(Tel, μ, E) .* (1 .- FermiDirac(Tel, μ, E))
end

function dFDdE!(tmp, Tel, μ, E)
    @. tmp = -1/(Constants.kB*Tel) * FermiDirac(Tel, μ, E) * (1 - FermiDirac(Tel, μ, E))
end

function dFDdE_mul!(tmp, Tel, μ, E)
    @. tmp = tmp * -1/(Constants.kB*Tel) * FermiDirac(Tel, μ, E) * (1 - FermiDirac(Tel, μ, E))
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
function dFDdT(Tel, μ, E)
    return (E-μ)/(Constants.kB*Tel^2) .* FermiDirac(Tel, μ, E) .* (1 .- FermiDirac(Tel, μ, E))
end

function dFDdT!(tmp, Tel, μ, E)
    @. tmp = (E-μ)/(Constants.kB*Tel^2) .* FermiDirac(Tel, μ, E) .* (1 .- FermiDirac(Tel, μ, E))
end

function dFDdT_mul!(tmp, Tel, μ, E)
    @. tmp = tmp * (E-μ)/(Constants.kB*Tel^2) .* FermiDirac(Tel, μ, E) .* (1 .- FermiDirac(Tel, μ, E))
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
function dFDdμ(Tel, μ, E)
    return 1/(Constants.kB*Tel) .* FermiDirac(Tel, μ, E) .* (1 .- FermiDirac(Tel, μ, E))
end

function dFDdμ!(tmp, Tel, μ, E)
    @. tmp = 1/(Constants.kB*Tel) * FermiDirac(Tel, μ, E) * (1 - FermiDirac(Tel, μ, E))
end

function dFDdμ_mul!(tmp, Tel, μ, E)
    @. tmp = tmp * 1/(Constants.kB*Tel) * FermiDirac(Tel, μ, E) * (1 - FermiDirac(Tel, μ, E))
end


function boltzmann_E_excitation(f, sim, E_mag, DOS)
    kgrid = sim.stucture.bandstructure.E_to_k(sim.structure.egrid)
    γ = bessel_gamma(E_mag, sim)
    κ = boltzmann_screening(f, kgrid, sim)
    fspl = get_interpolant(sim.structure.egrid, f)
    dfdt = zeros(length(f))
    for i in eachindex(f)
        k = kgrid[i]
        prefac = sim.electronicdistribution.Ω * pi / Constants.ħ / k
        for l in -3:3
            E1 = sim.structure.egrid[i] + l*sim.laser.hv
            k1 = sim.stucture.bandstructure.E_to_k(E1)
            dfrac = DOS(E1) / k1
            F = pauli_excitation_blocking(fspl, sim.structure.egrid[i], E1)
            mom_change = boltzmann_E_momentumintegral(l, k, k1, κ, γ, sim)
            dfdt[i] += dfrac * F * mom_change
        end
        dfdt[i] *= prefac
    end
    return dfdt 
end

function pauli_excitation_blocking(f,E,E_prime)
    return (f(E_prime) * (1-f(E))) - (f(E) * (1-f(E_prime)))
end

function boltzmann_E_momentumintegral(l, k, k1, κ, γ, sim)
    int(u,p) = u* excitation_matrix(sim, u, κ)* average_Bessel(l, γ*u) * boltzmann_step(u, k, k1)
    prob = IntegralProblem(int, 0.0, 2*k)
    sol = solve(prob, HCubatureJL(initdiv=100), abstol=1e-4, reltol=1e-4)
    return sol.u
end

function average_Bessel(order, value)
    int(u,p) = besselj(order, value*u)^2
    prob = IntegralProblem(int, -1, 1)
    sol = solve(prob, HCubatureJL(initdiv=50), abstol=1e-5, reltol=1e-5)
    return 1/2 * sol.u
end

function bessel_gamma(E_mag, sim)
    return Constants.q * E_mag / (sim.electronicdistribution.me * photon_energytofrequency(sim.laser.hv))
end

function boltzmann_step(q, k, k1)
    val = (q^2 + k^2 - k1^2) / (2*k*q)
    if -1 ≤ val ≤ 1
        return 1.0
    else 
        return 0.0
    end
end

"""
    ee_collision_integral(egrid, fgrid, N, v0; Constants.ħ=1.0)

Compute the electron-electron collision term

    (df/dt)_{e-e}

for every energy in `egrid` using sampled data on an energy grid.

The distribution `f` is assumed to be sampled on the same monotonically
increasing energy grid `egrid`. The density-of-states factor `N` is a spline
or other callable interpolation object. The delta function is used to
eliminate `xi'` analytically via `xi' = eps + epsp - xi`, and values of `f`
outside the grid are taken to be zero.

Arguments
- `egrid`: Energy grid shared by `fgrid`.
- `fgrid`: Occupation values sampled on `egrid`.
- `N`: Density-of-states spline or callable interpolation object.
- `v0`: Interaction strength.
- `Constants.ħ`: Reduced Planck constant.

Returns
- A vector containing `(df/dt)_{e-e}` evaluated at each point in `egrid`.
"""
function ee_collision_integral(sim, DOS, fgrid, v0)
    egrid = sim.structure.egrid
    out = similar(fgrid, promote_type(eltype(egrid), eltype(fgrid), Float64))
    cache = EECollisionIntegralCache(egrid)
    ee_collision_integral!(out, sim, fgrid, DOS, v0, cache)
    return out
end

"""
    ee_collision_integral!(out, egrid, fgrid, N, v0; Constants.ħ=1.0)

Fill `out` with the collision term evaluated at every energy in `egrid`.
"""
mutable struct RethfeldCollisionIntegralCache{T}
    ngrid::Vector{T}
    kgrid::Vector{T}
    inner::Vector{Vector{T}}
    epsp_values::Vector{Vector{T}}
    me_eff::Float64
end

mutable struct OnoCollisionIntegralCache{T}
    ngrid::Vector{T}
    kgrid::Vector{T}
    inner::Vector{Vector{T}}
    epsp_values::Vector{Vector{T}}
    v0::Float64
end

function EECollisionIntegralCache(egrid::AbstractVector, val, type::Symbol)
    T = promote_type(eltype(egrid), Float64)
    n = length(egrid)
    nslots = max(Threads.nthreads(), Threads.maxthreadid())
    if type == :Rethfeld
        return RethfeldCollisionIntegralCache{T}(zeros(T, n), zeros(T, n), 
        [zeros(T, n) for _ in 1:nslots], 
        [zeros(T, n) for _ in 1:nslots], 
        val)
    elseif type == :Ono
        return OnoCollisionIntegralCache{T}(zeros(T, n), zeros(T, n), 
        [zeros(T, n) for _ in 1:nslots], 
        [zeros(T, n) for _ in 1:nslots], 
        val)
    end
end

function ee_collision_integral!(out, sim, fgrid, DOS, val, type::Symbol)
    cache = EECollisionIntegralCache(sim.structure.egrid, val, type)
    return ee_collision_integral!(out, sim, fgrid, DOS, cache)
end

function ee_collision_integral!(out, sim, fgrid, DOS, cache::RethfeldCollisionIntegralCache)
    egrid = sim.structure.egrid
    n = length(egrid)
    
    sample_on_grid!(cache.ngrid, DOS, egrid)
    sample_on_grid!(cache.kgrid, sim.structure.bandstructure.E_to_k, egrid)

    prefactor = pi^3 ./ (Constants.ħ*cache.kgrid) 
    κ = boltzmann_screening(fgrid, sim, cache.kgrid, cache.me_eff)
    #@inbounds for k in eachindex(egrid)
    Threads.@threads for k in eachindex(egrid)
        tid = Threads.threadid()
        inner = cache.inner[tid]
        epsp_values = cache.epsp_values[tid]
        f_eps = fgrid[k]
        k_eps = cache.kgrid[k]

        @inbounds for i in eachindex(egrid)
            #epsp = egrid[i]
            f_epsp = fgrid[i]
            n_epsp = cache.ngrid[i]
            k_epsp = cache.kgrid[i]

            for j in eachindex(egrid)
                idx_xip = k + i - j

                if idx_xip < 1 || idx_xip > n
                    inner[j] = zero(eltype(inner))
                    continue
                end

                f_xi = fgrid[j]
                n_xi = cache.ngrid[j]
                k_xi = cache.kgrid[j]
                n_xip = cache.ngrid[idx_xip]
                f_xip = fgrid[idx_xip]
                k_xip = cache.kgrid[idx_xip]

                a = (k_eps^2 + k_epsp^2 - k_xi^2 - k_xip^2) / (2*k_xi*k_xip)
                b = (k_eps*k_epsp) / (k_xi*k_xip)
                if abs(a) > 1 + b
                    inner[j] = zero(eltype(inner))
                    continue
                end

                loss = f_eps * f_epsp * (1 - f_xi) * (1 - f_xip)
                gain = (1 - f_eps) * (1 - f_epsp) * f_xi * f_xip

                Δk_min = max(abs(k_eps - k_xi), abs(k_xip - k_epsp))
                Δk_max = min(k_eps + k_xi, k_xip + k_epsp)
                Mee2_integral = ee_matrix_element_integral(Δk_min, Δk_max, κ)

                inner[j] = (n_xi/k_xi) * (n_xip/k_xip) * (gain - loss) * Mee2_integral
            end

            epsp_values[i] = (n_epsp/k_epsp) * integration_algorithm(inner, egrid)
        end

        out[k] = prefactor[k] * integration_algorithm(epsp_values, egrid)
    end

    return out 
end

function ee_collision_integral!(out, sim, fgrid, DOS, cache::OnoCollisionIntegralCache)
    egrid = sim.structure.egrid
    n = length(egrid)

    prefactor = 2pi / Constants.ħ * cache.v0^2
    sample_on_grid!(cache.ngrid, DOS, egrid)
    #@inbounds for k in eachindex(egrid)
    Threads.@threads for k in eachindex(egrid)
        tid = Threads.threadid()
        inner = cache.inner[tid]
        epsp_values = cache.epsp_values[tid]
        #eps = egrid[k]
        f_eps = fgrid[k]
        #n_eps = cache.ngrid[k]

        @inbounds for i in eachindex(egrid)
            #epsp = egrid[i]
            f_epsp = fgrid[i]
            n_epsp = cache.ngrid[i]

            for j in eachindex(egrid)
                idx_xip = k + i - j
                if idx_xip < 1 || idx_xip > n
                    inner[j] = zero(eltype(inner))
                    continue
                end

                f_xi = fgrid[j]
                n_xi = cache.ngrid[j]
                n_xip = cache.ngrid[idx_xip]
                f_xip = fgrid[idx_xip]

                loss = f_eps * f_epsp * (1 - f_xi) * (1 - f_xip)
                gain = (1 - f_eps) * (1 - f_epsp) * f_xi * f_xip

                inner[j] = n_xi * n_xip * (gain - loss)
            end

            epsp_values[i] = n_epsp * integration_algorithm(inner, egrid)
        end

        out[k] = prefactor * integration_algorithm(epsp_values, egrid)
    end
    return out
end

function sample_on_grid!(dest::AbstractVector, f, grid::AbstractVector)
    @assert length(dest) == length(grid) "destination and grid must have the same length"
    @inbounds for i in eachindex(grid)
        dest[i] = f(grid[i])
    end
    return dest
end

# closed-form ∫_{Δk_min}^{Δk_max} |M_ee(Δk,κ)|² dΔk, avoiding numerical quadrature
# inside the O(n³) loop. |M_ee(Δk,κ)|² = (q²/ϵ0)² / (Δk²+κ²)²
function ee_matrix_element_integral(Δk_min, Δk_max, κ)
    F(x) = x/(2κ^2*(x^2+κ^2)) + atan(x/κ)/(2κ^3)
    return (Constants.q^2/Constants.ϵ0)^2 * (F(Δk_max) - F(Δk_min))
end

function boltzmann_screening(f, sim, kgrid, me_eff)
    prefac = Constants.q^2 * Constants.me * me_eff / (Constants.ϵ0*Constants.ħ^2*pi^2)
    f_spl = get_interpolant(kgrid, f)
    int(u,p) = f_spl(u)
    prob = IntegralProblem(int, minimum(kgrid), maximum(kgrid))
    sol = solve(prob, HCubatureJL(initdiv=100), abstol=1e-4, reltol=1e-4)
    return sqrt(prefac * sol.u)
end

function electron_electron_matrix(Δk, κ)
    frac1 = Constants.q^2 / Constants.ϵ0
    frac2 = 1/ (Δk^2 + κ^2)
    return (frac1*frac2)^2
end

function boltzmann_E_electronphonon()
    kgrid = sim.stucture.bandstructure.E_to_k(sim.structure.egrid)
    κ = boltzmann_screening(f, kgrid, sim)
    γ = bessel_gamma(E_mag, sim)
    fspl = get_interpolant(sim.structure.egrid, f)
    gspl = get_interpolant(sim.structure.egrid, g)
    dfdt = zeros(length(f))
    for i in eachindex(f)
        k = kgrid[i]
        prefac = 2 * sim.electronicdistribution.Ω * pi^3 / Constants.ħ / k
        int(u,p) = boltzmann_epscatter_int(i, sim, κ, DOS, fspl, gspl, γ)
        prob = IntegralProblem(int, 0.0, sim.phononicdistribution.ED)
        sol = solve(prob, HCubatureJL(initdiv=10), abstol=1e-3, reltol=1e-3)
        dfdt[i] = sol.u * prefac
    end
    return dfdt
end

function boltzmann_epscatter_int(Eq, sim, κ, DOS, f, g, γ)
    q = Eq / (Constants.ħ * sim.phononicdistribution.cs)
    M = electron_phonon_matrix(sim, Eq, q, κ)
    D = sim.phononicdistribution.DOS_ph(Eq) / q
    tmp = 0.0
    for l in -3:3
        for i in [-1, 1]
            E1 = E + (i*Eq) - l*sim.laser.hv
            k1 = sim.structure.bandstructure.E_to_k(E1)
            D1 = DOS(E1) / k1
            F = pauli_epscatter_blocking(f, g, i, E, E1)
            J = average_Bessel(l, γ*q)
            Ξ = boltzmann_step(q, k, k1)
            tmp += M * D * D1 * F * J * Ξ
        end
    end
    return tmp
end

function electron_phonon_matrix(sim, Eq, q, κ)
    frac1 = Constants.q^2 / (2*Constants.ϵ0 * sim.electronicdistribution.Ω)
    frac2 = Eq / (q^2 + κ^2)
    return frac1 * frac2
end

function pauli_epscatter_blocking(f, g, pm, E, E1)
    return f(E1) * (1 - f(E)) * (g(Eq) + 1/2 + (pm * 1/2)) - f(E) * (1 - f(E1)) * (g(Eq) + 1/2 - (pm * 1/2))
end

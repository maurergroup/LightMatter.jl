
using Roots

function find_chemicalpotential(N_target::Float64, Tel::Float64, DOS, egrid::Vector{Float64},
                                tmp, tmp2, μ_init::Float64 = 0.0)
    # n(μ) - N_target: the function to zero
    function residual(μ)
        return get_thermalparticles(tmp2, μ, Tel, DOS, egrid) - N_target
    end

    # Exact derivative dn/dμ — plugs in your existing function directly
    function dresidual(μ)
        return p_μ(tmp, μ, Tel, DOS, egrid)
    end

    # Newton with automatic fallback to bisection if it struggles
    return find_zero((residual, dresidual), μ_init, Order8(), atol=1e-5, rtol=1e-5)
end

"""
    get_thermalparticles(μ::Float64, Tel::Float64, DOS::spl, egrid::Vector{Float64})
    
    Calculates number of electrons assuming a Fermi Dirac distribution

    # Arguments
    - 'μ': Chemical potential
    - 'Tel': Temperature of the electronic bath
    - 'DOS': Density-of-states of the system
    - 'egrid': Energy grid the simulation is solved over

    # Returns
    - The number of electrons
"""
function get_thermalparticles(μ, Tel, DOS, egrid)
    f = DOS(egrid) .* FermiDirac(Tel,μ,egrid)
    return integration_algorithm(f, egrid)
end

function get_thermalparticles(tmp, μ, Tel, DOS, egrid)
    @. tmp = DOS(egrid)
    FermiDirac_mul!(tmp, Tel, ForwardDiff.value(μ), egrid)
    return integration_algorithm(tmp, egrid)
end
"""
    get_noparticles(Dis::Vector{Float64}, DOS::spl, egrid::Vector{Float64})
    
    Calculates number of electrons in the given distribution

    # Arguments
    - 'Dis': Electronic distribution
    - 'DOS': Density-of-states of the system
    - 'egrid': Energy grid the simulation is solved over

    # Returns
    - The number of electrons
"""
function get_noparticles(Dis, DOS, egrid)
    f = DOS(egrid) .* Dis
    return integration_algorithm(f,egrid)
end

function get_noparticles(tmp, Dis, DOS, egrid)
    @. tmp = DOS(egrid)
    tmp .*= Dis
    return integration_algorithm(tmp,egrid)
end
"""
    p_T(μ::Float64, Tel::Float64, DOS::spl, egrid::Vector{Float64})
    
    Calculates the change in the number of particles of a Fermi Dirac distribution with respect to temperature

    # Arguments
    - 'μ': Electronic distribution
    - 'Tel': The electronic bath temperature
    - 'DOS': Density-of-states of the system
    - 'egrid': Energy grid the simulation is solved over

    # Returns
    - The value of dn/dT
"""
function p_T(μ, Tel, DOS, egrid)
    f = dFDdT(Tel,μ,egrid) .* DOS(egrid)
    return integration_algorithm(f, egrid)
end

function p_T(tmp, μ, Tel, DOS, egrid)
    @. tmp = DOS(egrid)
    dFDdT_mul!(tmp, Tel, μ, egrid)
    return integration_algorithm(tmp, egrid)
end
"""
    p_μ(μ::Float64, Tel::Float64, DOS::spl, egrid::Vector{Float64})
    
    Calculates the change in the number of particles of a Fermi Dirac distribution with respect to chemical potential

    # Arguments
    - 'μ': Electronic distribution
    - 'Tel': The electronic bath temperature
    - 'DOS': Density-of-states of the system
    - 'egrid': Energy grid the simulation is solved over

    # Returns
    - The value of dn/dμ
"""
function p_μ(μ, Tel, DOS, egrid)
    f = dFDdμ(Tel,μ,egrid) .* DOS(egrid)
    return integration_algorithm(f, egrid)
end

function p_μ(tmp, μ, Tel, DOS, egrid)
    @. tmp = DOS(egrid)
    dFDdμ_mul!(tmp, Tel, μ, egrid)
    return integration_algorithm(tmp, egrid)
end
"""
    get_internalenergy(Dis::Vector{Float64}, DOS::spl, egrid::Vector{Float64})
    
    Calculates the internal energy of the given electronic distribution

    # Arguments
    - 'Dis': Electronic distribution
    - 'DOS': Density-of-states of the system
    - 'egrid': Energy grid the simulation is solved over

    # Returns
    - The internal energy of the given distribution
"""
function get_internalenergy(Dis, DOS, egrid)
    f = DOS(egrid) .* Dis .* egrid
    return integration_algorithm(f, egrid)
end

function get_internalenergy(tmp, Dis, DOS, egrid)
    @. tmp = DOS(egrid)
    tmp .*= Dis
    tmp .*= egrid
    return integration_algorithm(tmp, egrid)
end
"""
    c_T(μ::Float64, Tel::Float64, DOS::spl, egrid::Vector{Float64})
    
    Calculates the change in internal energy of a Fermi Dirac distribution with respect to temperature

    # Arguments
    - 'μ': Electronic distribution
    - 'Tel': The electronic bath temperature
    - 'DOS': Density-of-states of the system
    - 'egrid': Energy grid the simulation is solved over

    # Returns
    - The value of dU/dT
"""
function c_T(μ, Tel, DOS, egrid)
    f = dFDdT(Tel,μ,egrid) .* DOS(egrid) .* egrid
    return integration_algorithm(f, egrid)
end

function c_T(tmp, μ, Tel, DOS, egrid)
    @. tmp = DOS(egrid)
    dFDdT_mul!(tmp, Tel, μ, egrid)
    tmp .*= egrid
    return integration_algorithm(tmp, egrid)
end
"""
    c_μ(μ::Float64, Tel::Float64, DOS::spl, egrid::Vector{Float64})
    
    Calculates the change in internal energy of a Fermi Dirac distribution with respect to chemical potential

    # Arguments
    - 'μ': Electronic distribution
    - 'Tel': The electronic bath temperature
    - 'DOS': Density-of-states of the system
    - 'egrid': Energy grid the simulation is solved over

    # Returns
    - The value of dU/dμ
"""
function c_μ(μ, Tel, DOS, egrid)
    f = dFDdμ(Tel,μ,egrid) .* DOS(egrid) .* egrid
    return integration_algorithm(f, egrid)
end

function c_μ(tmp, μ, Tel, DOS, egrid)
    @. tmp = DOS(egrid)
    dFDdμ_mul!(tmp, Tel, μ, egrid)
    tmp .*= egrid
    return integration_algorithm(tmp, egrid)
end
"""
    integration_algorithm(y::AbstractVector, x::AbstractVector)
    
    Performs numerical integration using Simpson's rule

    # Arguments
    - 'y': Vector of function values to integrate
    - 'x': Vector of x values corresponding to y

    # Returns
    - The integral value
"""
function integration_algorithm(y::AbstractVector, x)
    h = x[2] - x[1]
    n = length(y)
    
    # Simple Simpson's rule for quick computation if we have enough points
    if n >= 3
        # Direct Simpson's 1/3 rule
        integral = zero(promote_type(eltype(y), typeof(h)))
        
        # Apply Simpson's 1/3 to all pairs of intervals
        if mod(n - 1, 2) == 0  # Even number of intervals
            @inbounds for i in 1:2:n-2
                integral += (h / 3) * (y[i] + 4*y[i+1] + y[i+2])
            end
            return integral
        else  # Odd number of intervals
            # Simpson's 1/3 on first part
            @inbounds for i in 1:2:n-3
                integral += (h / 3) * (y[i] + 4*y[i+1] + y[i+2])
            end
            # Trapezoidal on last
            integral += (h / 2) * (y[n-1] + y[n])
            return integral
        end
    end
    
    return 0.0
end

function get_particlenumber(DOS, egrid, offset)
    if DOS isa Vector{spl}
        pn = zeros(length(DOS))
        for i in eachindex(DOS)
            pn[i] = get_thermalparticles(offset[i],1e-16,DOS[i], egrid)
        end
        return pn
    else
        return get_thermalparticles(offset[1],1e-16, DOS, egrid)
    end
end

function get_particlenumber(DOS::Missing, egrid, offset)
    return missing
end
"""
    find_chemicalpotential_ORIGINAL(no_part::Float64, Tel::Float64, DOS::spl, egrid::Vector{Float64})
    
    Original implementation - kept for reference. Has allocation issues.
"""
function find_chemicalpotential_ORIGINAL(no_part, Tel, DOS, egrid, μ0)
    # PROBLEM: This creates a closure that captures variables
    # Each solve iteration may allocate through ForwardDiff
    f(u,p) = no_part - get_thermalparticles(u, Tel, DOS, egrid)
    return solve(NonlinearProblem(f, μ0);alg=SimpleNewtonRaphson(), abstol=1e-3, reltol=1e-3).u
end
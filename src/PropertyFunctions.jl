"""
    find_chemicalpotential(no_part::Float64, Tel::Float64, DOS::spl, egrid::Vector{Float64})
    
    Determines the chemical potential at the current temperature

    # Arguments
    - 'no_part': Float64 of particles in the thermal electronic system
    - 'Tel': Temperature of the electronic bath
    - 'DOS': Density-of-states of the system
    - 'egrid': Energy grid the simulation is solved over

    # Returns
    - The value of the chemical potential
"""
function find_chemicalpotential(no_part, Tel, DOS, egrid)
    f(u,p) = no_part - get_thermalparticles(u, Tel, DOS, egrid)
    return solve(NonlinearProblem(f, 0.0);alg=SimpleNewtonRaphson(), abstol=1e-8, reltol=1e-8).u
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
    # Use spline interpolation for the distribution for smoother integration
    f = DOS(egrid) .* Dis .* egrid
    return integration_algorithm(f, egrid)
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
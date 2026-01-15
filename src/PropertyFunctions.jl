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
    return solve(NonlinearProblem(f, 0.0);alg=SimpleKlement(), abstol=1e-4, reltol=1e-4).u
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
    f = i -> DOS(i) * FermiDirac(Tel,μ,i)
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
#= function get_noparticles(Dis, DOS, egrid)
    f = DOS(egrid) .* Dis
    return integration_algorithm(f,egrid)
end =#

function get_noparticles(Dis, DOS, egrid)
    # Use spline interpolation for the distribution for smoother integration
    Dis_spl = get_interpolant(egrid, Dis)
    f = i -> DOS(i) * Dis_spl(i)
    return integration_algorithm(f, egrid)
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
    f = i -> dFDdT(Tel,μ,i) * DOS(i)
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
    f = i -> dFDdμ(Tel,μ,i) * DOS(i)
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
    Dis_spl = get_interpolant(egrid, Dis)
    f = i -> DOS(i) * Dis_spl(i) * i
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

function integration_algorithm(y::Function, x) 
    return adaptive_high_order(y, x[1], x[end])
end

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
"""
    adaptive_high_order(f, a, b; tol=1e-14, maxdepth=100)

Adaptive integration with aggressive refinement for maximum accuracy.
Uses multiple convergence criteria and deeper recursion.
"""
function adaptive_high_order(f, a, b; tol=1e-14, maxdepth=100)
    total_evals = [0]  # Track function evaluations
    
    function integrate_recursive(a, b, depth, fa, fb)
        total_evals[1] += 13  # GK-15 uses 15 points, minus 2 endpoints
        
        result, error = gauss_kronrod_15(f, a, b)
        
        # Multiple stopping criteria for accuracy
        if depth >= maxdepth
            return result, error
        end
        
        # Stricter convergence check
        if error < tol * max(abs(result), 1e-12) && error < tol
            return result, error
        end
        
        # Split and recurse
        mid = (a + b) / 2
        fmid = f(mid)
        
        left_result, left_error = integrate_recursive(a, mid, depth + 1, fa, fmid)
        right_result, right_error = integrate_recursive(mid, b, depth + 1, fmid, fb)
        
        return left_result + right_result, sqrt(left_error^2 + right_error^2)
    end
    
    fa, fb = f(a), f(b)
    result, error = integrate_recursive(a, b, 0, fa, fb)
    return result#, error, total_evals[1]
end

"""
    gauss_kronrod_15(f, a, b)

15-point Gauss-Kronrod rule (7-point Gauss embedded).
Higher order than GK-21 for very smooth functions.
"""
function gauss_kronrod_15(f, a, b)
    # 7-point Gauss nodes
    xg = [0.0, 0.4058451513773972, 0.7415311855993945, 0.9491079123427585]
    wg = [0.4179591836734694, 0.3818300505051189, 0.2797053914892767, 0.1294849661688697]
    
    # 15-point Kronrod nodes and weights
    xk = [0.0000000000000000, 0.2077849550078985, 0.4058451513773972,
          0.5860872354676911, 0.7415311855993945, 0.8648644233597691,
          0.9491079123427585, 0.9914553711208126]
    wk = [0.2094821410847278, 0.2044329400752989, 0.1903505780647854,
          0.1690047266392679, 0.1406532597155259, 0.1047900103222502,
          0.0630920926299786, 0.0229353220105292]
    
    c = (b - a) / 2
    d = (a + b) / 2
    
    result_g = wg[1] * f(d)
    result_k = wk[1] * f(d)
    
    for i in 2:length(xk)
        fplus = f(d + c * xk[i])
        fminus = f(d - c * xk[i])
        result_k += wk[i] * (fplus + fminus)
        
        if i == 3 || i == 5 || i == 7  # Gauss points
            idx = (i == 3) ? 2 : (i == 5 ? 3 : 4)
            result_g += wg[idx] * (fplus + fminus)
        end
    end
    
    return c * result_k, abs(c * (result_k - result_g))
end
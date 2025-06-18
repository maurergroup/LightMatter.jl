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
#= function find_chemicalpotential(no_part, Tel, DOS, egrid)
    f(u,p) = no_part - get_thermalparticles(u, Tel, DOS, egrid)
    return solve(NonlinearProblem(f, 0.0); abstol=1e-3, reltol=1e-3).u
end =#

function find_chemicalpotential(no_part, Tel, DOS, egrid)::Float64
    #noe = ForwardDiff.value(no_part)
    temp = ForwardDiff.value(Tel)
    f(u,p) = no_part - get_thermalparticles(u, temp, DOS, egrid)
    return sol = solve(NonlinearProblem(f, 0.0); abstol=1e-3, reltol=1e-3).u
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
    return extended_Bode(DOS(egrid) .* FermiDirac(Tel,μ,egrid), egrid)
end

#= function get_thermalparticles(μ::ForwardDiff.Dual, Tel::SparseConnectivityTracer.GradientTracer, DOS, egrid)
    μ = ForwardDiff.value(μ)
    return extended_Bode(DOS(egrid) .* FermiDirac(Tel,μ,egrid), egrid)
end =#
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
    return extended_Bode(Dis.*DOS(egrid),egrid)
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
    return extended_Bode(dFDdT(Tel,μ,egrid) .* DOS(egrid), egrid)
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
    return extended_Bode(dFDdμ(Tel,μ,egrid) .* DOS(egrid), egrid)
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
    return extended_Bode(Dis .* DOS(egrid) .* egrid, egrid)
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
    return extended_Bode(dFDdT(Tel,μ,egrid) .* DOS(egrid) .* egrid, egrid)
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
    return extended_Bode(dFDdμ(Tel,μ,egrid) .* DOS(egrid) .* egrid, egrid)
end
"""
    extended_Bode(y::Vector{Float64}, x::Vector{Float64})
    
    Performs numerical integration on a grid using the higher order Bode's method.
    Will integrate from end to end of the x vector

    # Arguments
    - 'y': The spectrum on a grid to be integrated
    - 'x': The grid the spectrum is on w

    # Returns
    - The integration value across the range
"""
function extended_Bode(y, x)
    n = length(x)
    h = x[2]-x[1]  # The spacing between points
    integral = 0.0

    # Determine how many intervals to use for each rule
    integral_limit = n - 1

    # Apply Composite 3/8 Rule for the remaining intervals
    integral += (2h / 45) * (
        7 * y[1] + 7 * y[integral_limit+1] +
        32 * sum(y[2:2:integral_limit]) +
        12 * sum(y[3:4:integral_limit]) +
        14 * sum(y[4:3:integral_limit])
    )
    return integral
end
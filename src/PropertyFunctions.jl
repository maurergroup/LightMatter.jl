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
    return sol = solve(NonlinearProblem(f, 0.0), SimpleNewtonRaphson(); abstol=1e-12, reltol=1e-12).u
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
    f = i -> DOS(egrid[i]) * FermiDirac(Tel,μ,egrid[i])
    return Bode_rule(f, egrid)
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
    f = i -> Dis[i] * DOS(egrid[i])
    return Bode_rule(f,egrid)
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
    f = i -> dFDdT(Tel,μ,egrid[i]) * DOS(egrid[i])
    return Bode_rule(f, egrid)
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
    f = i -> dFDdμ(Tel,μ,egrid[i]) * DOS(egrid[i])
    return Bode_rule(f, egrid)
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
    f = i -> Dis[i] * DOS(egrid[i]) * egrid[i]
    return Bode_rule(f, egrid)
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
    f = i -> dFDdT(Tel,μ,egrid[i]) * DOS(egrid[i]) * egrid[i]
    return Bode_rule(f, egrid)
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
    f = i -> dFDdμ(Tel,μ,egrid[i]) * DOS(egrid[i]) * egrid[i]
    return Bode_rule(f, egrid)
end
"""
    Bode_rule(y::Vector{Float64}, x::Vector{Float64})
    
    Performs numerical integration on a grid using the higher order Boole's method.
    Will integrate from end to end of the x vector

    # Arguments
    - 'y': The spectrum on a grid to be integrated
    - 'x': The grid the spectrum is on w

    # Returns
    - The integration value across the range
"""
function Bode_rule(y::Vector{Float64}, x)
    N = length(x)
    n = (N-1) ÷   4
    h = x[2] - x[1]

    integral = 0.0
    @inbounds for i in 1:n
        idx = 4 * (i - 1) + 1
        integral += (2h / 45) * (
            7y[idx] +
            32y[idx + 1] +
            12y[idx + 2] +
            32y[idx + 3] +
            7y[idx + 4]
        )
    end

    return integral
end

function Bode_rule(y::Function, x)
    N = length(x)
    n = (N-1) ÷ 4
    h = x[2] - x[1]

    integral = 0.0
    @inbounds for i in 1:n
        idx = 4 * (i - 1) + 1
        integral += (2h / 45) * (
            7y(idx) +
            32y(idx + 1) +
            12y(idx + 2) +
            32y(idx + 3) +
            7y(idx + 4)
        )
    end

    return integral
end
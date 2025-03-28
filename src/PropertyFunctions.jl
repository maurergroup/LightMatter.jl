"""
    find_chemicalpotential(no_part::Float64,Tel::Float64,DOS::spl,egrid::Vector{<:Real})
    Sets up and solves the non-linear problem of determing the chemical potential at the current 
    electronic temperature.
"""
function find_chemicalpotential(no_part::Float64,Tel::Float64,DOS::spl,egrid::Vector{<:Real})
    f(u) = no_part - get_thermalparticles(u,Tel,DOS,egrid)
    return solve(ZeroProblem(f,0.0),Order1();atol=1e-3,rtol=1e-3)
end
"""
    get_thermalparticles(μ::Float64,Tel::Float64,DOS::spl,egrid::Vector{<:Real})
    Generates the number of thermal particles at a given temperature and chemical potential
"""
function get_thermalparticles(μ::Float64,Tel::Float64,DOS::spl,egrid::Vector{<:Real})
    return extended_Bode(DOS(egrid).*FermiDirac(Tel,μ,egrid),egrid)
end
"""
    get_noparticles(Dis::Vector{<:Real},DOS::spl,egrid::Vector{<:Real})
    Determines the number of particles in any system using a distribution that is on the
    same grid as the energy grid.
"""
function get_noparticles(Dis::Vector{<:Real},DOS::spl,egrid::Vector{<:Real})
    return extended_Bode(Dis.*DOS(egrid),egrid)
end
"""
    p_T(μ::Float64,Tel::Float64,DOS::spl,egrid::Vector{<:Real})
    Calculates the change in number of particles due to a change in temperature
"""
function p_T(μ::Float64,Tel::Float64,DOS::spl,egrid::Vector{<:Real})
    return extended_Bode(dFDdT(Tel,μ,egrid).*DOS(egrid),egrid)
end
"""
    p_μ(μ::Float64,Tel::Float64,DOS::spl,egrid::Vector{<:Real})
    Calculates the change in number of particles due to a change in chemical potential
"""
function p_μ(μ::Float64,Tel::Float64,DOS::spl,egrid::Vector{<:Real})
    return extended_Bode(dFDdμ(Tel,μ,egrid).*DOS(egrid),egrid)
end
"""
    get_internalenergy(Dis::Vector{<:Real},DOS::spl,egrid::Vector{<:Real})
    Determines the number of particles in any system using a distribution that is on the
    same grid as the energy grid.
"""
function get_internalenergy(Dis::Vector{<:Real},DOS::spl,egrid::Vector{<:Real})
    return extended_Bode(Dis.*DOS(egrid).*egrid,egrid)
end
"""
    c_T(μ::Float64,Tel::Float64,DOS::spl,egrid::Vector{<:Real})
    Calculates the change in internal energy due to a change in temperature
"""
function c_T(μ::Float64,Tel::Float64,DOS::spl,egrid::Vector{<:Real})
    return extended_Bode(dFDdT(Tel,μ,egrid).*DOS(egrid).*egrid,egrid)
end
"""
    c_μ(μ::Float64,Tel::Float64,DOS::spl,egrid::Vector{<:Real})
    Calculates the change in internal energy due to a change in chemical potential
"""
function c_μ(μ::Float64,Tel::Float64,DOS::spl,egrid::Vector{<:Real})
    return extended_Bode(dFDdμ(Tel,μ,egrid).*DOS(egrid).*egrid,egrid)
end
"""
    extended_Bode(y::Vector{<:Real}, x::Vector{<:Real})
    A high-order integration method beyond Simpson's rule. It is exact for polynomials up to fifth order
    and uses the quadrature rule for weighted evaluations. 
"""
function extended_Bode(y::Vector{<:Real}, x::Vector{<:Real})
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
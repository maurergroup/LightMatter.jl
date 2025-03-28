"""
    Generates an interpolation object with extrapolation from any two vectors of reals. The
    extrapolation returns the value of the boundaries. This should be suitable for DOS that are
    constant at the calculated boundaries and electronic distributions whose energy range is wide
    enough to capture all thermal and non-thermal behaviour.
"""
get_interpolate(xvals,yvals) = DataInterpolations.LinearInterpolation(yvals,xvals,extrapolation = ExtrapolationType.Constant)
"""
    Sets up and solves the non-linear problem of determing the chemical potential at the current 
    electronic temperature.
"""
function find_chemicalpotential(no_part::Float64,Tel::Float64,DOS::spl,kB::Float64,egrid)::Float64
    f(u) = no_part - get_thermalparticles(u,Tel,DOS,kB,egrid)
    return solve(ZeroProblem(f,0.0),Order1();atol=1e-3,rtol=1e-3)
end

function get_thermalparticles(μ::Float64,Tel::Float64,DOS::spl,kB::Float64,egrid)::Float64
    return extended_Bode(DOS(egrid).*FermiDirac(Tel,μ,kB,egrid),egrid)
end
"""
    Determines the number of particles in any system using an interpolation of the system and
    the DOS of the system.
"""
function get_noparticles(Dis::Vector{<:Real},DOS::spl,egrid)
    return extended_Bode(Dis.*DOS(egrid),egrid)
end

function p_T(μ::Float64,Tel::Float64,DOS::spl,kB::Float64,egrid)
    return extended_Bode(dFDdT(kB,Tel,μ,egrid).*DOS(egrid),egrid)
end

function p_μ(μ::Float64,Tel::Float64,DOS::spl,kB::Float64,egrid)
    return extended_Bode(dFDdμ(kB,Tel,μ,egrid).*DOS(egrid),egrid)
end
"""
    Determines the internal energy of any system using an interpolation of that system and the
    DOS of the system.
"""
function get_internalenergy(Dis::Vector{<:Real},DOS::spl,egrid)
    return extended_Bode(Dis.*DOS(egrid).*egrid,egrid)
end

function c_T(μ::Float64,Tel::Float64,DOS::spl,kB::Float64,egrid::Vector{<:Real})
    return extended_Bode(dFDdT(kB,Tel,μ,egrid).*DOS(egrid).*egrid,egrid)
end

function c_μ(μ::Float64,Tel::Float64,DOS::spl,kB::Float64,egrid)
    return extended_Bode(dFDdμ(kB,Tel,μ,egrid).*DOS(egrid).*egrid,egrid)
end

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
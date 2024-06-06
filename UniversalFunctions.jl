@kwdef struct Simulation_Misc <: simulation
    ERange::Vector{Float64}
    int_lb::Float64
    int_ub::Float64
end

struct GlobalConstants <: simulation
    kB::Float64
    ϵ0::Float64
    hbar::Float64
end
"""
    Returns the internal energy of the given function, requires an interpolation
    object of the function
"""
function get_internal_energy_int(func,mp::MaterialParameters,misc::Simulation_Misc)
    int(u,p) = func(u)*mp.DOS(u)*u
    prob=IntegralProblem(int,misc.int_lb,misc.int_ub)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-5,abstol=1e-5)
    return sol.u
end
"""
    Returns the internal energy of the given function from the Fermi-Dirac
    distriubtion
"""
function get_internal_energy(μ,Tel,mp::MaterialParameters,misc::Simulation_Misc,gc)
    int(u,p) = FermiDirac(u,μ,Tel,gc)*mp.DOS(u)*u
    prob=IntegralProblem(int,misc.int_lb,misc.int_ub)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-5,abstol=1e-5)
    return sol.u
end
"""
    Returns the number of electrons in the given function, requires an interpolation
    object of the function
"""
function get_noparticles_int(func,mp::MaterialParameters,misc::Simulation_Misc)
    int(u,p) = func(u)*mp.DOS(u)
    prob=IntegralProblem(int,misc.int_lb,misc.int_ub)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-5,abstol=1e-5)
    return sol.u
end
"""
    Returns the number of electrons in the given function, using the Fermi Dirac
    distribution
"""
function get_noparticles(μ,Tel,mp::MaterialParameters,misc::Simulation_Misc,gc)
    int(u,p) = FermiDirac(u,μ,Tel,gc)*mp.DOS(u)
    prob=IntegralProblem(int,misc.int_lb,misc.int_ub)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-5,abstol=1e-5)
    return sol.u
end
"""
    Returns an interpolation object of whatever data you sent to it with an
    extrapolation of 'nearest' meaning any calls outside the range return the boundary
    value
"""
function get_interpolate(xvals::Vector{Real},yvals::Vector{Real})
    return Spline1D(xvals,yvals,bc="nearest")
end
"""
    Generates an interpolation object that represents the DOS and shifted to the
    new Fermi energy
"""
function generate_DOS(File::String,FE,n)
    TotalDOS::Matrix{Float64}=readdlm(File,skipstart=3)
    return get_interpolate(TotalDOS[:,1].+FE,TotalDOS[:,2].*n)
end
"""
    The Fermi energy in this model is defined as the energy difference between the
    top and bottom of the valence band. This function finds that value assuming the
    top of the valence band is set to 0.
"""
function get_FermiEnergy(File)
    TotalDOS::Matrix{Float64}=readdlm(File,skipstart=3)
    nonzero=findfirst(TotalDOS[:,2] .== filter(x -> x!=0,TotalDOS[:,2])[1])
    return abs(TotalDOS[nonzero,1])
end

function find_chemicalpotential(no_part,Tel,μ,mp,misc)
    f(x) = no_part - get_noparticles(x,Tel,mp,misc)
    prob=ZeroProblem(f,μ)
    return solve(prob,Order16();atol=1e-5,reltol=1e-5)
end
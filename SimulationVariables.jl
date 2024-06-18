"""
    Determines the Fermi energy of the system. This is defined as the distance between the bottom and
    top of the valence band. This is used within Fermi Liquid Theory Relaxation Time where it is scaled
    based on the Fermi Energy^2. In all other places, the Fermi Energy is defined as 0.0.
"""
function get_FermiEnergy(File::String)
    TotalDOS::Matrix{Real}=readdlm(File,skipstart=3)
    Nonzero = findfirst(!=(0.0),TotalDOS[:,2])
    return abs(TotalDOS[Nonzero,1])
end
"""
    Converts a file location for the DOS into an interpolation object. It assumes that the DOS file
    is in units of states/atom and therefore scales the number of states by the number of atoms/nm(n).
"""
function generate_DOS(File::String,n)
    TotalDOS::Matrix{<:Real}=readdlm(File,skipstart=3)
    return get_interpolate(TotalDOS[:,1],TotalDOS[:,2].*n)
end
"""
    Generates an interpolation object with extrapolation from any two vectors of reals. The
    extrapolation returns the value of the boundaries. This should be suitable for DOS that are
    constant at the calculated boundaries and electronic distributions whose energy range is wide
    enough to capture all thermal and non-thermal behaviour.
"""
get_interpolate(xvals::Vector{Float64},yvals::Vector{Float64}) = Spline1D(xvals,yvals,bc="nearest")
@register_symbolic get_interpolate(xvals::Vector{Num},yvals::Vector{Num})::Spline1D
"""
    A callback function used to update the chemical potential with temperature. Is used when 
    Simulation.ParameterAPproximation.ChemicalPotential == true. ctx is a tuple holding the DOS 
     and number of particles.
"""
function update_chempot!(integ,u,p,ctx::Tuple{Spline1D,Real})
    integ.p[p.μ] = find_chemicalpotential(ctx[2],integ.u[u.Tel],integ.p[p.μ],ctx[1],integ.p[p.kB])
end
"""
    Sets up and solves the non-linear problem of determing the chemical potential at the current 
    electronic temperature.
"""
function find_chemicalpotential(no_part::Real,Tel::Real,μ::Real,DOS::Spline1D,kB::Real)
    f(u,p) = no_part - get_thermalparticles(u,Tel,DOS,kB)
    return solve(NonlinearProblem(f,μ),SimpleKlement();abstol=1e-3,reltol=1e-3).u
end
"""
    Finds the number of particles within a thermal system using the DOS and current temperature.
    This is solved to initially find no_part using T=1e-16 (as 0.0 leads to a discontinuity)
    and during finding the chemical potential. 
"""
function get_thermalparticles(μ::Real,Tel::Real,DOS::Spline1D,kB::Real)
    p=(μ,Tel,kB,DOS)
    int(u,p) = get_thermalparticles_int(u,p)
    return solve(IntegralProblem(int,(μ-10,μ+10),p),HCubatureJL(initdiv=10);abstol=1e-3,reltol=1e-3).u
end

function get_thermalparticles(μ::ForwardDiff.Dual,Tel::Real,DOS::Spline1D,kB::Real)
    p=(μ,Tel,kB,DOS)
    int(u,p) = get_thermalparticles_int(u,p)
    return solve(IntegralProblem(int,(ForwardDiff.value(μ)-10,ForwardDiff.value(μ)+10),p),HCubatureJL(initdiv=10);reltol=1e-3,abstol=1e-3).u
end
"""
    The integrand for finding the number of thermal particles. Using a parameter tuple with 
    components p=(μ,Tel,kB,DOS)
"""
function get_thermalparticles_int(u::Real,p::Tuple{Real,Real,Real,Spline1D})
    return FermiDirac(u,p[1],p[2],p[3])*p[4](u)
end
"""
    Determines the number of particles in any system using an interpolation of the system and
    the DOS of the system.
"""
function get_noparticles(μ::Real,Dis::Spline1D,DOS::Spline1D)
    int(u,p) = Dis(u) * DOS(u)
    return solve(IntegralProblem(int,(μ-10,μ+10),p),HCubatureJL(initdiv=2);reltol=1e-3,abstol=1e-3).u
end
"""
    Determines the internal energy of any system using an interpolation of that system and the
    DOS of the system.
"""
function get_internalenergy(μ::Real,Dis::Spline1D,DOS::Spline1D)
    int(u,p) = Dis(u) * DOS(u) * u
    return solve(IntegralProblem(int,(μ-10,μ+10),p),HCubatureJL(initdiv=2);reltol=1e-3,abstol=1e-3).u
end
@register_symbolic get_internalenergy(μ::Num,Dis::Spline1D,DOS::Spline1D)
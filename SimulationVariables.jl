"""
    Determines the Fermi energy of the system. This is defined as the distance between the bottom and
    top of the valence band. This is used within Fermi Liquid Theory Relaxation Time where it is scaled
    based on the Fermi Energy^2. In all other places, the Fermi Energy is defined as 0.0.
"""
function get_FermiEnergy(File::String)
    TotalDOS::Matrix{Float64}=readdlm(File)
    Nonzero = findfirst(!=(0.0),TotalDOS[:,2])
    return abs(TotalDOS[Nonzero,1])
end
"""
    Converts a file location for the DOS into an interpolation object. It assumes that the DOS file
    is in units of states/atom and therefore scales the number of states by the number of atoms/nm(n).
"""
function generate_DOS(File::String,n)
    TotalDOS::Matrix{Float64}=readdlm(File)
    return get_interpolate(TotalDOS[:,1],TotalDOS[:,2].*n)
end
"""
    Generates an interpolation object with extrapolation from any two vectors of reals. The
    extrapolation returns the value of the boundaries. This should be suitable for DOS that are
    constant at the calculated boundaries and electronic distributions whose energy range is wide
    enough to capture all thermal and non-thermal behaviour.
"""
get_interpolate(xvals,yvals) = DataInterpolations.LinearInterpolation(yvals,xvals,extrapolate=true)
"""
    Sets up and solves the non-linear problem of determing the chemical potential at the current 
    electronic temperature.
"""
function find_chemicalpotential(no_part::Float64,Tel::Float64,DOS::spl,kB::Float64,FE::Float64,n0::Float64)::Float64
    f(u) = no_part - get_thermalparticles(u,Tel,DOS,kB,FE,n0)
    return solve(ZeroProblem(f,0.0),Order1();atol=1e-3,rtol=1e-3)
end

function get_thermalparticles(μ::Float64,Tel::Float64,DOS::spl,kB::Float64,FE::Float64,n0::Float64)::Float64
    int_neg(u,p) = DOS(u)*(1/(exp((u-μ)/(kB*Tel))+1)-1)
    int_pos(u,p) = DOS(u)/(exp((u-μ)/(kB*Tel))+1)
    prob=IntegralProblem(int_neg,(-FE,0.0))
    prob2=IntegralProblem(int_pos,(0.0,20.0))
    uroot_neg = solve(prob,HCubatureJL(initdiv=10);abstol=1e-5,reltol=1e-5).u
    uroot_pos = solve(prob2,HCubatureJL(initdiv=10);abstol=1e-5,reltol=1e-5).u
    return uroot_neg+uroot_pos+n0
end
"""
    Determines the number of particles in any system using an interpolation of the system and
    the DOS of the system.
"""
function get_noparticlesspl(Dis::spl,DOS::spl,n0,FE)
    int_neg(u,p) = (Dis(u).-1)*DOS(u)
    int_pos(u,p) = Dis(u)*DOS(u)
    prob=IntegralProblem(int_neg,(-FE,0.0))
    prob2=IntegralProblem(int_pos,(0.0,20.0))
    uroot_neg = solve(prob,HCubatureJL(initdiv=10);abstol=1e-5,reltol=1e-5).u
    uroot_pos = solve(prob2,HCubatureJL(initdiv=10);abstol=1e-5,reltol=1e-5).u
    return uroot_neg+uroot_pos+n0
end

function get_n0(DOS,μ,FE)
    int(u,p) = DOS(u)
    prob=IntegralProblem(int,(-FE,μ))
    return solve(prob,HCubatureJL(initdiv=10),reltol=1e-5,abstol=1e-5).u
end

function p_T(μ::Float64,Tel::Float64,DOS::spl,kB::Float64)
    int(u,p) = dFDdT(kB,Tel,μ,u)*DOS(u)
    prob=IntegralProblem(int,(μ-(60*Tel/10000),μ+(60*Tel/10000)))
    return solve(prob,HCubatureJL(initdiv=10);abstol=1e-5,reltol=1e-5).u
end

function p_μ(μ::Float64,Tel::Float64,DOS::spl,kB::Float64)
    int(u,p) = dFDdμ(kB,Tel,μ,u)*DOS(u)
    prob=IntegralProblem(int,(μ-(60*Tel/10000),μ+(60*Tel/10000)))
    return solve(prob,HCubatureJL(initdiv=10);abstol=1e-5,reltol=1e-5).u
end
"""
    Determines the internal energy of any system using an interpolation of that system and the
    DOS of the system.
"""
function get_internalenergyspl(Dis::spl,DOS::spl,u0::Float64,FE)
    int_neg(u,p) = (Dis(u)-1)*DOS(u)*u
    int_pos(u,p) = Dis(u)*DOS(u)*u
    prob=IntegralProblem(int_neg,(-FE,0.0))
    prob2=IntegralProblem(int_pos,(0.0,20.0))
    u_neg = solve(prob,HCubatureJL(initdiv=10),reltol=1e-5,abstol=1e-5).u
    u_pos = solve(prob2,HCubatureJL(initdiv=10),reltol=1e-5,abstol=1e-5).u
    return u_neg+u_pos+u0
end

function get_u0(DOS,μ,FE)
    int(u,p) = DOS(u)*u
    prob=IntegralProblem(int,(-FE,μ))
    return solve(prob,HCubatureJL(initdiv=10),reltol=1e-5,abstol=1e-5).u
end

function c_T(μ::Float64,Tel::Float64,DOS::spl,kB::Float64)
    int(u,p) = dFDdT(kB,Tel,μ,u).*DOS(u)*u
    prob = IntegralProblem(int,(μ-(60*Tel/10000),μ+(60*Tel/10000)))
    return solve(prob,HCubatureJL(initdiv=10);reltol=1e-5,abstol=1e-5).u
end

function c_μ(μ::Float64,Tel::Float64,DOS::spl,kB::Float64)
    int(u,p) = dFDdμ(kB,Tel,μ,u).*DOS(u)*u
    prob = IntegralProblem(int,(μ-(60*Tel/10000),μ+(60*Tel/10000)))
    return solve(prob,HCubatureJL(initdiv=10);reltol=1e-5,abstol=1e-5).u
end
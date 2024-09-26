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
get_interpolate(xvals,yvals) = Spline1D(xvals,yvals,bc="nearest")
"""
    A callback function used to update the chemical potential with temperature. Is used when 
    Simulation.ParameterAPproximation.ChemicalPotential == true. ctx is a tuple holding the DOS 
     and number of particles.
"""
function update_chempotTTM!(integ,u,p,ctx::Tuple{Spline1D,Real})
    integ.p[p.μ] = find_chemicalpotential(ctx[2],integ.u[u.Tel],integ.p[p.μ],ctx[1],integ.p[p.kB])
end

function update_chempotAthEM!(integ,u,p,ctx::Tuple{Spline1D})
    integ.p[p.μ] = find_chemicalpotential(integ.u[u.n],integ.u[u.Tel],integ.p[p.μ],ctx[1],integ.p[p.kB])
end
"""
    Sets up and solves the non-linear problem of determing the chemical potential at the current 
    electronic temperature.
"""
function find_chemicalpotential(no_part::Real,Tel::Real,μ::Real,DOS::Spline1D,kB::Real)
    f(u) = no_part - get_thermalparticles(u,Tel,DOS,kB)
    return solve(ZeroProblem(f,μ),Order1();atol=1e-3,rtol=1e-3)
end

function get_thermalparticles(μ,Tel::Real,DOS::Spline1D,kB::Real)
    int(u,p) = FermiDirac(u,μ,Tel,kB)*DOS(u)
    return solve(IntegralProblem(int,(μ-20,μ+20)),HCubatureJL(initdiv=50);abstol=1e-6,reltol=1e-6).u
end
"""
    Determines the number of particles in any system using an interpolation of the system and
    the DOS of the system.
"""
function get_noparticlesspl(μ::Real,Dis::Spline1D,DOS::Spline1D,n0)
    int_neg(u,p) = (Dis(u).-1)*DOS(u)
    uroot_neg = solve(IntegralProblem(int_neg,(-Inf,0.0)),CubatureJLh();reltol=1e-3,abstol=1e-3).u
    int_pos(u,p) = Dis(u)*DOS(u)
    uroot_pos = solve(IntegralProblem(int_pos,(0.0,μ+10)),CubatureJLh();reltol=1e-3,abstol=1e-3).u
    return uroot_neg+uroot_pos+n0
end

function get_noparticles(Dis,DOS::Spline1D,egrid)
    integrand = Dis.*DOS(egrid)
    prob = SampledIntegralProblem(integrand,egrid)
    return solve(prob,SimpsonsRule()).u
end
@register_symbolic get_noparticles(Dis::AbstractVector,DOS::Spline1D,egrid::AbstractVector)

function get_n0(DOS,μ)
    int(u,p) = DOS(u)
    prob=IntegralProblem(int,(-Inf,μ))
    return solve(prob,HCubatureJL(initdiv=100),reltol=1e-5,abstol=1e-5).u
end

function p_T(μ::Real,Tel::Real,DOS::Spline1D,kB::Real)
    p = (μ,Tel,kB,DOS)
    int = BatchIntegralFunction(p_T_int,zeros(0))
    return solve(IntegralProblem(int,(μ-(6*Tel/10000),μ+(6*Tel/10000)),p),CubatureJLh();reltol=1e-5,abstol=1e-5).u
end
@register_symbolic p_T(μ::Num,Tel::Num,DOS::Spline1D,kB::Num)
function p_T_int(y,u,p)
    Threads.@threads for i in 1:length(u)
        @inbounds y[i] = dFDdT(p[3],p[2],p[1],u[i]).*p[4].(u[i])
    end
end

function p_μ(μ::Real,Tel::Real,DOS::Spline1D,kB::Real)
    p = (μ,Tel,kB,DOS)
    int = BatchIntegralFunction(p_μ_int,zeros(0))
    return solve(IntegralProblem(int,(μ-(6*Tel/10000),μ+(6*Tel/10000)),p),CubatureJLh();reltol=1e-5,abstol=1e-5).u
end
@register_symbolic p_μ(μ::Num,Tel::Num,DOS::Spline1D,kB::Num)
function p_μ_int(y,u,p)
    Threads.@threads for i in 1:length(u)
        @inbounds y[i] = dFDdμ(p[3],p[2],p[1],u[i]).*p[4].(u[i])
    end
end

"""
    Determines the internal energy of any system using an interpolation of that system and the
    DOS of the system.
"""
function get_internalenergyspl(μ::Real,Dis::Spline1D,DOS::Spline1D,u0::Real)
    int_neg(u,p) = (Dis(u).-1)*DOS(u)*u
    uroot_neg = solve(IntegralProblem(int_neg,(-Inf,0.0)),CubatureJLh();reltol=1e-3,abstol=1e-3).u
    int_pos(u,p) = Dis(u)*DOS(u)*u
    uroot_pos = solve(IntegralProblem(int_pos,(0.0,μ+10)),CubatureJLh();reltol=1e-3,abstol=1e-3).u
    return uroot_neg+uroot_pos+u0
end

function get_internalenergy(μ::Real,Tel::Real,DOS::Spline1D,kB::Real)
    p = (μ,Tel,kB,DOS)
    int = BatchIntegralFunction(internalenergy_int,zeros(0))
    return solve(IntegralProblem(int,(-Inf,10.0),p),CubatureJLh();reltol=1e-5,abstol=1e-5).u
end
@register_symbolic get_internalenergy(μ::Num,Tel::Num,DOS::Spline1D,kB::Num)
function internalenergy_int(y,u,p)
    Threads.@threads for i in 1:length(u)
        @inbounds y[i] = FermiDirac(p[2],p[1],p[3],u[i]).*p[4].(u[i]) *u[i]
    end
end

function get_internalenergy_grid(Dis,DOS,egrid)
    integrand = Dis.*DOS(egrid).*egrid
    prob = SampledIntegralProblem(integrand,egrid)
    return solve(prob,SimpsonsRule()).u
end
@register_symbolic get_internalenergy_grid(Dis::AbstractVector,DOS::Spline1D,egrid::AbstractVector)
function get_u0(DOS,μ)
    int(u,p) = DOS(u)*u
    prob=IntegralProblem(int,(-Inf,μ))
    return solve(prob,HCubatureJL(initdiv=100),reltol=1e-5,abstol=1e-5).u
end

function c_T(μ::Real,Tel::Real,DOS::Spline1D,kB::Real)
    p = (μ,Tel,kB,DOS)
    int = BatchIntegralFunction(c_T_int,zeros(0))
    return solve(IntegralProblem(int,(μ-(6*Tel/10000),μ+(6*Tel/10000)),p),CubatureJLh();reltol=1e-5,abstol=1e-5).u
end
@register_symbolic c_T(μ::Num,Tel::Num,DOS::Spline1D,kB::Num)
function c_T_int(y,u,p)
    Threads.@threads for i in 1:length(u)
        @inbounds y[i] = dFDdT(p[3],p[2],p[1],u[i]).*p[4].(u[i]) *u[i]
    end
end

function c_μ(μ::Real,Tel::Real,DOS::Spline1D,kB::Real)
    p = (μ,Tel,kB,DOS)
    int = BatchIntegralFunction(c_μ_int,zeros(0))
    return solve(IntegralProblem(int,(μ-(6*Tel/10000),μ+(6*Tel/10000)),p),CubatureJLh();reltol=1e-5,abstol=1e-5).u
end
@register_symbolic c_μ(μ::Num,Tel::Num,DOS::Spline1D,kB::Num)
function c_μ_int(y,u,p)
    Threads.@threads for i in 1:length(u)
        @inbounds y[i] = dFDdμ(p[3],p[2],p[1],u[i]).*p[4].(u[i]) *u[i]
    end
end
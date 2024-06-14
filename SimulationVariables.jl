function generate_DOS(File::String,n)
    TotalDOS::Matrix{Float64}=readdlm(File,skipstart=3)
    return get_interpolate(TotalDOS[:,1],TotalDOS[:,2].*n)
end

get_interpolate(xvals::Vector{Float64},yvals::Vector{Float64}) = Spline1D(xvals,yvals,bc="nearest")

function update_chempot!(integ,u,p,ctx::Tuple{Spline1D,Float64})
    integ.p[p.μ] = find_chemicalpotential(ctx[2],integ.u[u.Tel],integ.p[p.μ],ctx[1],integ.p[p.kB])
end

function find_chemicalpotential(no_part::Float64,Tel::Float64,μ::Float64,DOS::Spline1D,kB::Float64)
    f(u,p) = no_part - get_thermalparticles(u,Tel,DOS,kB)
    return solve(NonlinearProblem(f,μ),SimpleKlement();abstol=1e-3,reltol=1e-3).u
end

#= function get_thermalparticles_int(y::Vector{<:Real},u::Vector{<:Real},p::Tuple{Real,Float64,Float64,Spline1D})
    n=Threads.nthreads()
    Threads.@threads for i in 1:n
        @inbounds y[i:n:end] .= FermiDirac.(@view(u[i:n:end]),p[1],p[2],p[3]).*p[4].(@view(u[i:n:end]))
    end
end

function get_thermalparticles(μ::Float64,Tel::Float64,DOS::Spline1D,kB::Float64)
    p=(μ,Tel,kB,DOS)
    int = BatchIntegralFunction(get_thermalparticles_int,zeros(0))
    return solve(IntegralProblem(int,(μ-10,μ+10),p),QuadGKJL();abstol=1e-3).u
end

function get_thermalparticles(μ::ForwardDiff.Dual,Tel::Float64,DOS::Spline1D,kB::Float64)
    p=(μ,Tel,kB,DOS)
    int = BatchIntegralFunction(get_thermalparticles_int,Array{ForwardDiff.Dual}(undef,0))
    return solve(IntegralProblem(int,(ForwardDiff.value(μ)-10,ForwardDiff.value(μ)+10),p),QuadGKJL();reltol=1e-3,abstol=1e-3).u
end =#

function get_thermalparticles_int(u::Real,p::Tuple{Real,Float64,Float64,Spline1D})
    return FermiDirac(u,p[1],p[2],p[3])*p[4](u)
end

function get_thermalparticles(μ::Float64,Tel::Float64,DOS::Spline1D,kB::Float64)
    p=(μ,Tel,kB,DOS)
    int(u,p) = get_thermalparticles_int(u,p)
    return solve(IntegralProblem(int,(μ-10,μ+10),p),HCubatureJL(initdiv=10);abstol=1e-3,reltol=1e-3).u
end

function get_thermalparticles(μ::ForwardDiff.Dual,Tel::Float64,DOS::Spline1D,kB::Float64)
    p=(μ,Tel,kB,DOS)
    int(u,p) = get_thermalparticles_int(u,p)
    return solve(IntegralProblem(int,(ForwardDiff.value(μ)-10,ForwardDiff.value(μ)+10),p),HCubatureJL(initdiv=10);reltol=1e-3,abstol=1e-3).u
end

function get_noparticles(μ::Float64,Dis::Spline1D,DOS::Spline1D)
    int(u,p) = Dis(u) * DOS(u)
    return solve(IntegralProblem(int,(μ-10,μ+10),p),HCubatureJL(initdiv=2);reltol=1e-3,abstol=1e-3).u
end

function get_internalenergy(μ::Float64,Dis::Spline1D,DOS::Spline1D)
    int(u,p) = Dis(u) * DOS(u) * u
    return solve(IntegralProblem(int,(μ-10,μ+10),p),HCubatureJL(initdiv=2);reltol=1e-3,abstol=1e-3).u
end

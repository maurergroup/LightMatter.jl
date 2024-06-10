function generate_DOS(File::String,FE,n)
    TotalDOS::Matrix{Float64}=readdlm(File,skipstart=3)
    return get_interpolate(TotalDOS[:,1].+FE,TotalDOS[:,2].*n)
end

get_interpolate(xvals::Vector{Float64},yvals::Vector{Float64}) = Spline1D(xvals,yvals,bc="nearest")

function find_chemicalpotential(no_part::Float64,Tel::Float64,μ::Float64,DOS::Spline1D,kB::Float64)
    f(u,p) = no_part - get_noparticles(u,Tel,DOS,kB)
    sol = solve(NonlinearProblem(f,μ),SimpleKlement();abstol=1e-3,reltol=1e-3)
    return sol.u
end

function get_noparticles_int(y::Vector{Float64},u::Vector{Float64},p::Tuple{Float64,Float64,Float64,Spline1D})
    n=Threads.nthreads()
    Threads.@threads for i in 1:n
        @inbounds y[i:n:end] .= FermiDirac.(@view(u[i:n:end]),p[1],p[2],p[3]).*p[4](@view(u[i:n:end]))
    end
end

function get_noparticles_int(y::Vector{ForwardDiff.Dual},u::Vector{Float64},p::Tuple{ForwardDiff.Dual,Float64,Float64,Spline1D})
    n=Threads.nthreads()
    Threads.@threads for i in 1:n
        @inbounds y[i:n:end] .= FermiDirac.(@view(u[i:n:end]),p[1],p[2],p[3]).*p[4](@view(u[i:n:end]))
    end
end

function get_noparticles(μ::Float64,Tel::Float64,DOS::Spline1D,kB::Float64)
    p=(μ,Tel,kB,DOS)
    int = BatchIntegralFunction(get_noparticles_int,zeros(0))
    return solve(IntegralProblem(int,(μ-10,μ+10),p),QuadGKJL();abstol=1e-3).u
end

function get_noparticles(μ::ForwardDiff.Dual,Tel::Float64,DOS::Spline1D,kB::Float64)
    p=(μ,Tel,kB,DOS)
    int = BatchIntegralFunction(get_noparticles_int,Array{ForwardDiff.Dual}(undef,0))
    return solve(IntegralProblem(int,(ForwardDiff.value(μ)-10,ForwardDiff.value(μ)+10),p),QuadGKJL();reltol=1e-3,abstol=1e-3).u
end


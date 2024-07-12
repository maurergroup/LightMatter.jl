using BenchmarkTools,NonlinearSolve,Roots,Integrals,Dierckx,DelimitedFiles,StaticArrays,Plots

function generate_DOS(File::String,n)
    TotalDOS::Matrix{<:Real}=readdlm(File,skipstart=3)
    return get_interpolate(TotalDOS[:,1],TotalDOS[:,2].*n)
end

function find_relaxeddistribution(goal::Real,no_part::Real,DOS::Spline1D,kB::Real)
    f(u) = goal - find_temperatureandμ(u,no_part,DOS,kB)
    Temp = solve(ZeroProblem(f,1000.0),Order16();abstol=1e-5,reltol=1e-5)
    μ = find_chemicalpotential(no_part,Temp,0.0,DOS,kB)
    return Temp,μ
end

function find_temperatureandμ(Tel::Real,no_part::Real,DOS::Spline1D,kB::Real)
    μ = find_chemicalpotential(no_part,Tel,0.0,DOS,kB)
    return get_internalenergy(μ,Tel,DOS,kB)
end

function find_chemicalpotential(no_part::Real,Tel::Real,μ::Real,DOS::Spline1D,kB::Real)
    f(u) = no_part - get_thermalparticles(u,Tel,DOS,kB)
    return solve(ZeroProblem(f,μ),Order16();atol=1e-2,rtol=1e-2)
end

function get_thermalparticles(μ::Real,Tel::Real,DOS::Spline1D,kB::Real)
    p=(μ,Tel,kB,DOS)
    int = BatchIntegralFunction(get_thermalparticlesint,zeros(0))
    return solve(IntegralProblem(int,(μ-10,μ+10),p),CubatureJLh();abstol=1e-3,reltol=1e-3).u
end

function get_thermalparticlesint(y,u,p::Tuple{Real,Real,Real,Spline1D})
    Threads.@threads for i in 1:length(u)
        @inbounds y[i] = FermiDirac(p[2],p[1],p[3],u[i]).*p[4].(u[i])
    end
end

function get_internalenergy(μ::Real,Tel::Real,DOS::Spline1D,kB::Real)
    p = (μ,Tel,kB,DOS)
    int = BatchIntegralFunction(internalenergy_int,zeros(0))
    return solve(IntegralProblem(int,(μ-10,μ+10),p),CubatureJLh();reltol=1e-6,abstol=1e-6).u
end

function internalenergy_int(y,u,p)
    Threads.@threads for i in 1:length(u)
        @inbounds y[i] = FermiDirac(p[2],p[1],p[3],u[i]).*p[4].(u[i]) *u[i]
    end
end

FermiDirac(Tel,μ,kB,E) = 1 ./(exp.((E.-μ)./(kB*Tel)).+1)

function main()
    DOS = generate_DOS("DOS/Au_DOS.dat",59)
    no_part = get_thermalparticles(0.0,1e-18,DOS,8.617e-5)

    Tel=collect(range(500.0,5000.0,step=250.0))
    μ = find_chemicalpotential.(no_part,Tel,0.0,Ref(DOS),8.617e-5)
    goal = get_internalenergy.(μ,Tel,Ref(DOS),8.617e-5)
    CalcTel = zeros(length(Tel))
    Calcμ = zeros(length(Tel))
    for i in eachindex(Tel)
        println(Tel[i])
        CalcTel[i],Calcμ[i]=find_relaxeddistribution(goal[i],no_part,DOS,8.617e-5)
    end
    p1=plot(Tel,[Tel,CalcTel])
    p2=plot(μ,[μ,Calcμ])
    plot(p1,p2)
    #= Tel = 2500.0
    μ = find_chemicalpotential(no_part,Tel,0.0,DOS,8.617e-5)
    goal = get_internalenergy(μ,Tel,DOS,8.617e-5)
    println(Tel," ",μ)
    @btime find_relaxeddistribution($goal,$no_part,$DOS,8.617e-5) =#
end

using DelimitedFiles,Plots,Integrals,Interpolations

function tau_0(plasma)
    τ=128/(sqrt(3)*pi^2*plasma)
    return τ
end

function tau_ee(τ,μ,T_el,E,kB)
    τee=τ*μ^2 /((E-μ)^2 +(pi*kB*T_el)^2)
    return τee
end

function normalise(FDext,μ)
    int(u,p)=FDext(u)
    prob=IntegralProblem(int,μ-8,μ+8)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-5,abstol=1e-5)
    return sol.u
end

function convint(μ,E,FDext,σ)
    int(u,p)=gaussian(σ,E-u)*FDext(u)
    prob=IntegralProblem(int,μ-15,μ+15)
    sol=solve(prob,HCubatureJL(initdiv=10000);reltol=1e-10,abstol=1e-10)
    return sol.u
end

function gaussian(σ,ϵ)
    return 1/(σ*sqrt(2*pi))*exp(-ϵ^2/(2σ^2))
end

function gaus_nearestneighbour(FDext,E,σ,ω)
    temp=zeros(length(ω))
    for (i,W) in enumerate(ω)
        temp[i]=FDext(E+W)*gaussian(σ,W)
    end
    return sum(temp./sum(gaussian.(σ,ω)))
end

function lin_nearestneighbour(FDext,E,σ,ω)
    temp=zeros(length(ω))
    for (i,W) in enumerate(ω)
        temp[i]=FDext(E+W)*σ[i]
    end
    return sum(temp./sum(σ))
end

function storage()
    for (i,E) in enumerate(ERange)
        FD[i]=FDext(E)
        FDp1[i]=(FDext(E)+FDext(E+0.1)+FDext(E-0.1))/3
        FDp2[i]=(FDext(E)+FDext(E+0.2)+FDext(E-0.2))/3
        FDp5[i]=(FDext(E)+FDext(E+0.5)+FDext(E-0.5))/3
    end

    #ERange=range(FE-6,FE+6,step=0.001)
    FD=zeros(length(ERange))
    FDp1=zeros(length(ERange))
    FDp2=zeros(length(ERange))
    FDp3=zeros(length(ERange))
    FDp4=zeros(length(ERange))

    Threads.@threads for i in eachindex(ERange)
        FD[i]=FDext(ERange[i])
        FDp1[i]=lin_nearestneighbour(FDext,ERange[i],[1,3,1],collect(-0.1:0.05:0.05))
        FDp2[i]=lin_nearestneighbour(FDext,ERange[i],[1/3,1,3,1,1/3],collect(-0.1:0.05:0.1))
        FDp3[i]=lin_nearestneighbour(FDext,ERange[i],[1,3,1],collect(-0.1:0.1:0.1))
        FDp4[i]=lin_nearestneighbour(FDext,ERange[i],[1/3,1,3,1,1/3],collect(-0.2:0.1:0.2))
        #FDconv[i]=convint(FE,ERange[i],FDext,0.05)
    end
    #plot(ERange,[FD,FDconv])
    plot(ERange,[FD,FDp1,FDp2,FDp3,FDp4],title="Linear Averaging",
    label=["FD" "p1" "p2" "p3" "p4"])


    FD=FDext.(ERange)
    FDt1=lin_nearestneighbour.(Ref(FDext),ERange,Ref([1,3,1]),Ref(collect(-0.1:0.1:0.1)))
    FDt1int=Interpolations.interpolate(ERange,FDt1,SteffenMonotonicInterpolation())
    FDt1ext = extrapolate(FDt1int, Flat())
    FDt2=lin_nearestneighbour.(Ref(FDt1ext),ERange,Ref([1,3,1]),Ref(collect(-0.1:0.1:0.1)))
    FDt2int=Interpolations.interpolate(ERange,FDt2,SteffenMonotonicInterpolation())
    FDt2ext = extrapolate(FDt2int, Flat())
    FDt3=lin_nearestneighbour.(Ref(FDt2ext),ERange,Ref([1,3,1]),Ref(collect(-0.1:0.1:0.1)))
    FDt3int=Interpolations.interpolate(ERange,FDt3,SteffenMonotonicInterpolation())
    FDt3ext = extrapolate(FDt3int, Flat())

    plot(ERange,[FD,FDt1,FDt2,FDt3])
end

function DensityOfStates(File,FE,n)
    TotalDOS::Matrix{Float64}=readdlm(File,skipstart=3)
    DOS=Interpolations.interpolate(TotalDOS[:,1].+FE, TotalDOS[:,2]*n, SteffenMonotonicInterpolation())
    DOSext=extrapolate(DOS,Flat())
    return DOSext
end

function NumberOfElectrons(FDint,lb,ub,DOS)::Float64
    int(u,p)=DOS(u)*FDint(u)
    prob=IntegralProblem(int,lb,ub)
    sol=solve(prob,HCubatureJL(initdiv=100);reltol=1e-5,abstol=1e-5)
    return sol.u
end

function innerint(E,Eprime,DOS,FE)
    int(u,p)=DOS(u)*DOS(E+u-Eprime)
    prob=IntegralProblem(int,FE-E-Eprime,FE)
    sol=solve(prob,HCubatureJL(initdiv=2);abstol=1e-5,reltol=1e-5)
    return sol.u
end

function Aeschlimanntau(E,DOS,FE)
    int(u,p)=DOS(u)*innerint(E,u,DOS,FE)
    if E < FE
        prob=IntegralProblem(int,E,FE)
    else 
        prob=IntegralProblem(int,FE,E)
    end
    sol=solve(prob,HCubatureJL(initdiv=2);abstol=1e-4,reltol=1e-4)
    return sol.u
end

function matrixelements(DOS,plasma,FE,hbar,E)
    return sqrt(3)*pi/128/DOS(E)^2*hbar*plasma/FE^2
end

function NoEScaling(FDint,DOS,lb,ub)
    NoE=NumberOfElectrons(FDint,lb,ub,DOS)
    #TotNoE=NumberOfElectrons(neqext,FE-10,FE+10,DOS)
    return NoE
end

function DOSint(DOS,lb,ub,FDint)
    int(u,p)=DOS(u)*(1-FDint(u))
    prob=IntegralProblem(int,lb,ub)
    sol=solve(prob,HCubatureJL(initdiv=100);reltol=1e-5,abstol=1e-5)
    return sol.u
end

function main()
    IOFD=readdlm("NoneqTest.csv",skipstart=1)
    IOEnergy=IOFD[:,1].+9.9
    IODis=IOFD[:,2]
    neqint=Interpolations.interpolate(IOEnergy,abs.(IODis),SteffenMonotonicInterpolation())
    neqext = extrapolate(neqint, Flat())
    TotalFD=readdlm("../PhD Code/Gold/FermiDis/0fs_3.1_10F.csv",skipstart=1)
    TotEnergy=parse.(Float64,chop.(TotalFD[:,1],head=1,tail=1)).+9.9
    TotDis=parse.(Float64,chop.(TotalFD[:,2],tail=1))
    Totint=Interpolations.interpolate(TotEnergy,abs.(TotDis),SteffenMonotonicInterpolation())
    Totext=extrapolate(Totint, Flat())
    τ=tau_0(2.1835)
    T_el=3000
    kB=8.617e-5
    FE=9.9
    hbar=6.582e-1
    n=59
    plasma=5.8
    DOS=DensityOfStates("DOS/Au_DOS.dat",FE,n)
    ERange=range(FE-3,FE+3,step=0.01)
    NoEFD=zeros(length(ERange))
    NoDOS=zeros(length(ERange))
    tau=tau_ee.(τ,FE,T_el,ERange,kB)

    σ=0.1
    Threads.@threads for i in eachindex(ERange)
        NoEFD[i]=NoEScaling(neqext,DOS,ERange[i]-σ,ERange[i]+σ)
        NoDOS[i]=DOSint(DOS,ERange[i]-σ,ERange[i]+σ,neqext)
    end

    return ERange, tau.*NoDOS./NoEFD
end
ERange,Test = main()
plot(ERange,Test)
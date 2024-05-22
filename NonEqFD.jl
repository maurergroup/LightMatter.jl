using Integrals,Plots,SpecialFunctions,DelimitedFiles,Interpolations,Roots,LaTeXStrings
using OrdinaryDiffEq,BenchmarkTools

include("Structs.jl")
using .Structs

function tau_0(mp,sv)
    sv.τ=128/(sqrt(3)*pi^2*mp.ω)
end

function tau_ee(E,sv,gc)
    return sv.τ*sv.μ^2 /((E-sv.μ)^2 +(pi*gc.kB*sv.Tel[1])^2)
end

function DensityOfStates(mp,sv,Scale::Float64)
    TotalDOS=readdlm(mp.DOSFile,skipstart=3)
    DOS=Interpolations.interpolate(TotalDOS[:,1].+mp.FE, TotalDOS[:,2]*mp.n*Scale, SteffenMonotonicInterpolation())
    sv.DOS=DOS
    sv.lb=TotalDOS[1,1]+mp.FE
    sv.ub=TotalDOS[end,1]+mp.FE
end

function NumberOfElectrons(Tel,μ,sv,gc)::Float64
    int(u,p)=sv.DOS(u)*FermiDirac(u,Tel,μ,gc)
    prob=IntegralProblem(int,sv.lb,sv.ub)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-5,abstol=1e-5)
    return sol.u
end

function NumberOfElectrons(sv,gc,μ)::Float64
    int(u,p)=sv.DOS(u)*FermiDirac(u,sv.Tel[1],μ,gc)
    prob=IntegralProblem(int,sv.lb,sv.ub)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-5,abstol=1e-5)
    return sol.u
end

function ChemicalPotential(sv,gc)::Float64
    x0=sv.μ
    f(x)=sv.NumElec-NumberOfElectrons(sv,gc,x)
    prob=ZeroProblem(f,x0)
    sv.μ=solve(prob,Order16();atol=1e-3,rtol=1e-3)
end

function ChemicalPotential(Tel::Float64,μ::Float64,sv,gc)::Float64
    x0=μ
    f(x)=sv.NumElec-NumberOfElectrons(Tel,x,sv,gc)
    prob=ZeroProblem(f,x0)
    sv.μ=solve(prob,Order16();atol=1e-3,rtol=1e-3)
end

function FermiDirac(E::Float64,Tel,μ,gc)::Float64
    return 1/(exp((E-μ)/(gc.kB*Tel))+1)
end

function FermiDirac(E::Float64,sv,gc)::Float64
    return 1/(exp((E-sv.μ)/(gc.kB*sv.Tel[1]))+1)
end

function probint(sv,E_h::Float64,E_e::Float64,gc)::Float64
    return sv.DOS(E_h)*FermiDirac(E_h,sv,gc)*(1-FermiDirac(E_e,sv,gc))
end

function particlenumber(vec::Vector{Float64},sv,ERange::StepRangeLen)::Float64
    return sum(vec.*sv.DOS.(ERange).*(ERange[2]-ERange[1]))
end

function particleconstant(sv,l,gc)::Float64
    ERange=range(-4+sv.μ,4+sv.μ,step=0.01)
    elec=probint.(Ref(sv),ERange.-l.hv,ERange,Ref(gc))
    hol=probint.(Ref(sv),ERange,ERange.+l.hv,Ref(gc))
    f(x)=(particlenumber(elec*x,sv,ERange))-particlenumber(hol,sv,ERange)
    x0=1
    prob=ZeroProblem(f,x0)
    return solve(prob,Order8();atol=1e-3,rtol=1e-3)
end

function FDChange(E::Float64,sv,particle::Float64,gc,l)::Float64
    elec=probint(sv,E-l.hv,E,gc)
    hol=probint(sv,E,E+l.hv,gc)
    return (particle*elec).-hol
end

function integrand(l,t,τ)
    A=-4*log(2)/l.FWHM^2
    return exp(A*t^2)*exp(t*((1/τ)-(2*A*l.delay)))
end

function integral(l,τ,t)
    int(u,p)=integrand(l,t,τ)
    prob=IntegralProblem(int,0.0,t)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-10,abstol=1e-10)
    return sol.u
end

function erf1(A,B,t)
    return ((2*A*t)+B)/(2*sqrt(complex(A)))
end

function erf2(A,B)
    return B/(2*sqrt(complex(A)))
end

function exponential(A,B)
    return exp(-B^2/(4*A))
end

function solvedintegral(τ,t)
    A=-4*log(2)/l.FWHM^2
    B=(1/τ)-(2*A*l.delay)
    c1=exponential(A,B)
    c2=erf1(A,B,t)
    c3=erf2(A,B)
    return (sqrt(pi)*c1*(erfi(c2)-erfi(c3)))/(2*sqrt(complex(A)))
end

function prefactor(sv,gc,mp,l,E::Float64,particle::Float64,internalFD::Float64
    ,τ::Float64,t::Float64,B::Float64,C::Float64)

    num=4*FDChange(E,sv,particle,gc,l)*sqrt(log(2))*l.ϕ*mp.α*exp((B*l.delay^2)-(C^2/(4*B)))
    denom=l.FWHM*internalFD*2*sqrt(complex(B))

    return num/denom
end

function gaussian(σ,ϵ)
    return 1/(σ*sqrt(2*pi))*exp(-ϵ^2/(2σ^2))
end

function convint(μ,E,FDext,σ)
    int(u,p)=gaussian(σ,E-u)*FDext(u)
    prob=IntegralProblem(int,μ-15,μ+15)
    sol=solve(prob,HCubatureJL(initdiv=10000);reltol=1e-10,abstol=1e-10)
    return sol.u
end

function differential(sv,gc,l,mp,E::Float64,t::Float64,particle::Float64,τ::Float64,internalFD,B::Float64)
    C=(1/τ)+(8*l.delay*log(2)/l.FWHM^2)
    A=prefactor(sv,gc,mp,l,E,particle,internalFD,τ,t,B,C)
    f1=2*A*sqrt(complex(B))*exp(((2*B*t)+C)^2/(4*B)-(t/τ))/sqrt(pi)
    f2=A*exp(-t/τ)*(erfi(C/(2*sqrt(complex(B))))-erfi(((2*B*t)+C)/(2*sqrt(complex(B)))))/τ
    return f1+f2
end

function difftest(sv,gc,l,mp,Ttraj::Vector{Float64},t::Float64
    ,ERange::StepRangeLen,FD::Array{Float64},tstep::Float64)

    sv.Tel[1]=Ttraj[trunc(Int,t*10)+1]
    sv.μ=ChemicalPotential(sv,gc)
    particle=particleconstant(sv,l,gc)
    τ=tau_ee.(ERange,Ref(sv),Ref(gc))
    B=-4*log(2)/l.FWHM^2

    int(u,p)=FDChange(u,sv,particle,gc,l)*sv.DOS(u)*u
    prob=IntegralProblem(int,sv.μ-(3*l.hv),sv.μ+(3*l.hv))
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-10,abstol=1e-10)
    internalFD=sol.u

    Threads.@threads for i in eachindex(ERange)
        FD[i]+=differential(sv,gc,l,mp,ERange[i],t,particle,τ[i],internalFD,B)*tstep
        #= if ERange[i]<=sv.μ
            τ=excholerelax(ERange[i],sv,gc)
        elseif ERange[i]>sv.μ
            τ=excelecrelax(ERange[i],sv,gc)
        end
        dFDdt=differential(sv,gc,l,mp,ERange[i],t,particle,τ,internalFD,B)*tstep
        if isnan(dFDdt)
            return 0.0
        else
            FD[i]+=dFDdt
        end =#
    end

end

function Fermivector(ERange::StepRangeLen,sv,gc)
    Fermi=zeros(length(ERange))
    sv.μ=ChemicalPotential(sv,gc)
    Threads.@threads for E in eachindex(ERange)
        Fermi[E]=FermiDirac(ERange[E],sv,gc)
    end
    return Fermi
end

function main()
    l,mp,sim,gc,sv=Structs.parameterbuilder("InputFiles/Au_Input.txt")
    Ttraj=parse.(Float64,chop.(readdlm("../PhD Code/Gold/0D/0D_3.1_10F_5fsNoep.csv",skipstart=1)[:,2],tail=1))
    DensityOfStates(mp,sv,1.0)
    ERange=range(-4+mp.FE,4+mp.FE,step=0.001)
    
    sv.NumElec=NumberOfElectrons(0.0,mp.FE,sv,gc)
    sv.τ=tau_0(mp,sv)

    FD=zeros(length(ERange))
    tstop=300.0
    tstep=1.0
    Diffpoints=[350]
    Diffout=Array{Number}(undef,length(FD),length(Diffpoints))
    trange=range(0.0,tstop,step=tstep)
    counter=1
    for t in trange
        println(t)
        difftest(sv,gc,l,mp,Ttraj,t,ERange,FD,tstep)
        #= if t in Diffpoints
            Diffout[:,counter]=FD
            counter+=1
        end =#
    end
    #ThermFermi=Fermivector(ERange,sv,gc)
    #plot(ERange,Diffout,label=["300" "350" "400" "500"])
    #p2=plot(ERange,FD.+ThermFermi)
    #plot(p1,p2)
    output=hcat(ERange.-sv.μ,FD)
    writedlm("NoneqTest.csv",' ')
    io=open("NoneqTest.csv","a")
    for row in eachrow(output)
        println(io,join(row," "))
    end
    close(io)
end

main()
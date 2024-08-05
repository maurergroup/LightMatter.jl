using Integrals,DelimitedFiles,Dierckx,LaTeXStrings,JLD2,RecursiveArrayTools
using ForwardDiff,Roots,OrdinaryDiffEq,Unitful,NonlinearSolve

function DensityOfStates(File)
    TotalDOS::Matrix{Float64}=readdlm(File,skipstart=3)
    return Spline1D(TotalDOS[:,1].+FE,TotalDOS[:,2]*n,bc="nearest")
end

function get_internalenergy(μ::Real,Dis::Spline1D,DOS::Spline1D)
    int(u,p) = Dis(u) * DOS(u) * u
    sol = solve(IntegralProblem(int,(μ-10,μ+10)),HCubatureJL(initdiv=2);reltol=1e-3,abstol=1e-3).u
    return sol
end

function get_internalenergy(μ::Real,Tel::Real,DOS::Spline1D)
    int(u,p) = FermiDirac(u,μ,Tel) * DOS(u) * u
    sol = solve(IntegralProblem(int,(μ-10,μ+10)),HCubatureJL(initdiv=2);reltol=1e-3,abstol=1e-3).u
    return sol
end

function get_noparticles(μ::Real,Dis::Spline1D,DOS::Spline1D)
    int(u,p) = Dis(u) * DOS(u)
    return solve(IntegralProblem(int,(μ-10,μ+10)),HCubatureJL(initdiv=2);reltol=1e-3,abstol=1e-3).u
end

function get_FermiEnergy(File)
    TotalDOS::Matrix{Real}=readdlm(File,skipstart=3)
    Nonzero = findfirst(!=(0.0),TotalDOS[:,2])
    return abs(TotalDOS[Nonzero,1])
end

function FermiDirac(E::Float64,μ::Float64,Tel)::Float64
    return 1/(exp((E-μ)/(kB*Tel))+1)
end

function NumberOfElectrons(μ::Float64,Tel::Float64,DOS)::Float64
    int(u,p)=DOS(u)*FermiDirac(u,μ,Tel)
    prob=IntegralProblem(int,0.0,μ+10)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-10,abstol=1e-10)
    return sol.u
end

function ChemicalPotential(Tel,DOS,nopart)
    x0=FE
    f(x)=nopart-NumberOfElectrons(x,Tel,DOS)
    prob=ZeroProblem(f,x0)
    return solve(prob,Order16();atol=1e-10,rtol=1e-10)
end

function Cph(Tph) 
    uplim=θ./Tph #Generates an upperlimit for each z position in the Tph
    int(u,p)=(u^4*exp(u))/((exp(u)-1)^2)
    prob = IntegralProblem(int,0.0,uplim)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-10,abstol=1e-10) #sovles the integral within Cph for each position in Tph
    cph::Float64=(9*n*kB)*(Tph/θ)^3 *sol.u # Finds Cph for each z coordinate by broadcasting equation over Tph and integral
    return cph
end 

function tau_ee(τ,μ,Tel,E)
    τee=τ*μ^2 /((E-μ)^2 +(pi*kB*Tel)^2)
    return τee
end

function tau_th(g,Tel,Tph)
    τth=(γ*(Tel+Tph)/(2*g))
    return τth
end

function dFDdE(Tel::Float64,μ::Float64,E::Float64)::Real
    Numer=-exp((E-μ)/(kB*Tel))
    Denom=kB*Tel*(exp((E-μ)/(kB*Tel))+1)^2
    return Numer/Denom
end

function dFDdT(Tel::Float64,μ::Float64,E::Float64)::Real
    Numer=(E-μ)*exp((E-μ)/(kB*Tel))
    Denom=kB*Tel^2*(exp((E-μ)/(kB*Tel))+1)^2
    return Numer/Denom
end

function dFDdμ(Tel::Float64,μ::Float64,E::Float64)::Real
    Numer=exp((E-μ)/(kB*Tel))
    Denom=kB*Tel*(exp((E-μ)/(kB*Tel))+1)^2
    return Numer/Denom
end

function gep_nonlinear(Tel,μ,DOS)
    int(u,p)=DOS(u)^2*-dFDdE(Tel,μ,u)
    prob=IntegralProblem(int,(μ-(6*Tel/10000),μ+(6*Tel/10000)))
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-10,abstol=1e-10)
    λ=23e-6
    gep::Float64=sol.u*(pi*kB*λ)/(DOS(μ)*hbar)
    return gep
end

function GTelTph(Tel,Tph,μ,DOS)
    return gep_nonlinear(Tel,μ,DOS)*(Tel-Tph)
end

function athem_neq(ftot,DOS,tprime,eescatter,μ)
    excite = athem_excitation(ftot,DOS,tprime,μ)
    return excite .+ eescatter
end

function LaserStrength(tprime)
    return α*sqrt(4*log(2)/pi)*(ϕ/FWHM)*exp(-4*log(2)*((tprime-Delay)/FWHM)^2)
end

function athem_excitation(ftot,DOS,tprime,μ)
    Δfneqh = fgr_hole_generation(DOS,ftot)
    Δfneqe = fgr_electron_generation(DOS,ftot)
    pc_sf = fgr_particleconservation(DOS,Δfneqh,Δfneqe,μ)
    δ = LaserStrength(tprime)/fgr_excitation_internalenergy((Δfneqe*pc_sf) - Δfneqh,DOS,μ)
    return δ*((Δfneqe*pc_sf) - Δfneqh)
end

function fgr_hole_generation(DOS::Spline1D,ftotspl)
    return DOS.(egrid.+hv).*ftotspl.(egrid).*(1 .-ftotspl.(egrid.+hv))
end

function fgr_electron_generation(DOS::Spline1D,ftotspl)
    return DOS.(egrid.-hv).*ftotspl.(egrid.-hv).*(1 .-ftotspl.(egrid))
end

function fgr_particleconservation(DOS::Spline1D,fneqh::Vector{Float64},fneqe::Vector{Float64},μ::Real)
    elDis = Spline1D(egrid,fneqe,bc="nearest")
    hDis = Spline1D(egrid,fneqh,bc="nearest")
    f(u,p) = get_noparticles(μ,hDis,DOS) - u*get_noparticles(μ,elDis,DOS)
    return solve(NonlinearProblem(f,1.0),Klement();abstol=1e-3,reltol=1e-3).u
end

function fgr_excitation_internalenergy(Δfneq::Vector{<:Real},DOS::Spline1D,μ::Real)
    fneqspl = Spline1D(egrid,Δfneq,bc="nearest")
    return get_internalenergy(μ,fneqspl,DOS)
end

function athem_electronelectronscattering(ftot,Tel,μ,DOS,fneq,feq,τ,nopart)
    τee = tau_ee.(τ,μ,Tel,egrid)
    frel = find_relaxeddistribution(get_internalenergy(μ,ftot,DOS),nopart,DOS)
    return (fneq.+frel.-feq)./τee
end

function find_relaxeddistribution(goal::Real,no_part::Real,DOS::Spline1D)
    f(u,p) = goal - find_temperatureandμ(u,no_part,DOS)
    Temp = solve(NonlinearProblem(f,1000.0),Klement();abstol=1e-5,reltol=1e-5).u
    μ = ChemicalPotential(Temp,DOS,no_part)
    return FermiDirac.(egrid,μ,Temp)
end

function find_temperatureandμ(Tel::Real,no_part::Real,DOS::Spline1D)
    μ = ChemicalPotential(Tel,DOS,no_part)
    return get_internalenergy(μ,Tel,DOS)
end

function neqelectron_electrontransfer(eescatter,DOS,μ)
    disspl = Spline1D(egrid,eescatter,bc="nearest") 
    return get_internalenergy(μ,disspl,DOS)
end

function neqelectron_electronparticlechange(eescatter,DOS,μ)
    disspl = Spline1D(egrid,eescatter,bc="nearest") 
    return get_noparticles(μ,disspl,DOS)
end

function dUdT(DOS,Tel,μ)
    int(u,p) = DOS(u)*u*dFDdT(Tel,μ,u)
    prob = IntegralProblem(int,(μ-(6*Tel/10000),μ+(6*Tel/10000)))
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-10,abstol=1e-10)
    return sol.u
end

function dUdμ(DOS,Tel,μ)
    int(u,p) = DOS(u)*u*dFDdμ(Tel,μ,u)
    prob = IntegralProblem(int,(μ-(6*Tel/10000),μ+(6*Tel/10000)))
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-10,abstol=1e-10)
    return sol.u
end

function dNdT(DOS,Tel,μ)
    int(u,p) = DOS(u)*dFDdT(Tel,μ,u)
    prob = IntegralProblem(int,(μ-(6*Tel/10000),μ+(6*Tel/10000)))
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-10,abstol=1e-10)
    return sol.u
end

function dNdμ(DOS,Tel,μ)
    int(u,p) = DOS(u)*dFDdμ(Tel,μ,u)
    prob = IntegralProblem(int,(μ-(6*Tel/10000),μ+(6*Tel/10000)))
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-10,abstol=1e-10)
    return sol.u
end

function dTeldt(Tel,Tph,μ,DOS,eescatter,δn)
    G=GTelTph(Tel,Tph,μ,DOS)
    uee=neqelectron_electrontransfer(eescatter,DOS,μ)
    U=-G+uee
    cT = dUdT(DOS,Tel,μ)
    cμ = dUdμ(DOS,Tel,μ)
    pμ = dNdT(DOS,Tel,μ)
    pT = dNdμ(DOS,Tel,μ)
    return 1/(cT*pμ - pT*cμ)*(pμ*U-cμ*δn)
end

function Trajectory!(du,u,p,t) #currently no p arrays indexing
    println(t)
    p[1]=ChemicalPotential(u.x[2][1],p[3],u.x[1][1])
    feq = FermiDirac.(egrid,p[1],u.x[2][1])
    ftot = Spline1D(egrid,feq.+u.x[3][:],bc="nearest")
    eescatter = athem_electronelectronscattering(ftot,u.x[2][1],p[1],p[3],u.x[3][:],feq,p[2],u.x[1][1])

    du.x[1][1]=neqelectron_electronparticlechange(eescatter,p[3],p[1])
    du.x[2][1]=dTeldt(u.x[2][1],p[4],p[1],p[3],eescatter,du.x[1][1])
    du.x[3][:]= athem_neq(ftot,p[3],t,eescatter,p[1])
end

function Run0D(Input::Dict)
    τ = 0.546
    DOS=DensityOfStates(Input["DOS"])
    nopart=NumberOfElectrons(FE,0.0,DOS)
    μ=ChemicalPotential(300.0,DOS,nopart)
    tspan=(0.0,Input["SimEnd"])
    u0 = ArrayPartition([nopart],[300.0],zeros(length(egrid)))
    p=[μ,τ,DOS,300.0]
    prob=ODEProblem(Trajectory!,u0,tspan,p)
    println("Begin Solving")
    sol=solve(prob,AutoVern7(Rodas5(autodiff=false));reltol=1e-4,abstol=1e-4,saveat=0.1)
    @save Input["Output"] sol
end

function Units(Input)
    Input["g"]=ustrip(uconvert(u"eV/nm^3/fs/K", Input["g"]u"W/m^3/K"))
    Input["gamma"]=ustrip(uconvert(u"eV/nm^3/K^2",Input["gamma"]u"J/m^3/K^2"))
    Input["AbsLen"]=ustrip(uconvert(u"nm",Input["AbsLen"]u"m"))
    Input["AtomDens"]=ustrip(uconvert(u"nm^-3",Input["AtomDens"]u"m^-3"))
    Input["FWHM"]=ustrip(uconvert(u"fs",Input["FWHM"]u"s"))
    Input["Fluence"]=ustrip(uconvert(u"eV/nm^2",Input["Fluence"]u"J/m^2"))
    Input["LaserOff"]=ustrip(uconvert(u"fs",Input["LaserOff"]u"s"))
    Input["RTKappa"]=ustrip(uconvert(u"eV/nm/K/fs",Input["RTKappa"]u"W/m/K"))
    Input["SimEnd"]=ustrip(uconvert(u"fs",Input["SimEnd"]u"s"))
    Input["Length"]=ustrip(uconvert(u"nm",Input["Length"]u"m"))
    Input["dz"]=ustrip(uconvert(u"nm",Input["dz"]u"m"))
    Input["Plasma"]=Input["Plasma"]/(6.582e-1*2*pi)
    return Input
end

IO = readdlm("Au_Input.txt",'=',comments=true,comment_char='#')
for (x,y) in enumerate(IO[:,2])
    if typeof(y) == SubString{String}
        IO[x,2]=strip(IO[x,2],' ')
    end
end
Input=Dict(IO[i,1]=>IO[i,2] for i in 1:size(IO,1))
Input=Units(Input)
Input["AbsLen"]=1/Input["AbsLen"]
global const hv=Input["PhotEn"]
global const γ=Input["gamma"]
global const ϕ=Input["Fluence"]
global const α=Input["AbsLen"]
global const κrt=Input["RTKappa"]
global const FWHM=Input["FWHM"]
global const Delay=Input["LaserOff"]
global const n=Input["AtomDens"]
global const θ=Input["DebyeT"]
global const L=Input["Length"]
global const ne=Input["Freeelec"]
global const effmass=Input["Effmass"]
global const FE=get_FermiEnergy(Input["DOS"])
global const ω=Input["Plasma"]
global const egrid=range(FE-2*hv,FE+2*hv,step=0.01)
global const kB=8.617333e-5 #Boltzmann constant in eV K^-1
global const elecmass=1 #mass of electron in amu
global const eleccharge=1 #charge of electron in amu
global const ϵ0=8.065279e-9 #covnerted to amu dielectric constant of vacuum
global const hbar=6.582e-1 #eV fs
Run0D(Input) 

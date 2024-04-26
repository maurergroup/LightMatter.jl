module Equations

using DelimitedFiles,Roots,Integrals,Interpolations,ForwardDiff

function DensityOfStates(mp,sv)
    TotalDOS::Matrix{Float64}=readdlm(mp.DOSFile,skipstart=3)
    DOS=Interpolations.interpolate(TotalDOS[:,1].+mp.FE, TotalDOS[:,2]*mp.n, SteffenMonotonicInterpolation())
    sv.DOS=DOS
    sv.lb=TotalDOS[1,1]+mp.FE
    sv.ub=TotalDOS[end,1]+mp.FE
end

function DensityOfStates(mp,sv,Scale::Float64)
    TotalDOS::Matrix{Float64}=readdlm(mp.DOSFile,skipstart=4)
    DOS=Interpolations.interpolate(TotalDOS[:,1].+mp.FE, TotalDOS[:,2]*mp.n*Scale, SteffenMonotonicInterpolation())
    sv.DOS=DOS
    sv.lb=TotalDOS[1,1]+mp.FE
    sv.ub=TotalDOS[end,1]+mp.FE
end

function FermiDirac(E::Float64,Tel,μ,gc)::Float64
    return 1/(exp((E-μ)/(gc.kB*Tel))+1)
end

function FermiDirac(E::Float64,sv,gc)::Float64
    return 1/(exp((E-sv.μ)/(gc.kB*sv.Tel))+1)
end

function NumberOfElectrons(Tel,μ,sv,gc)::Float64
    int(u,p)=sv.DOS(u)*FermiDirac(u,Tel,μ,gc)
    prob=IntegralProblem(int,sv.lb,sv.ub)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-5,abstol=1e-5)
    return sol.u
end

function NumberOfElectrons(sv,gc,μ)::Float64
    int(u,p)=sv.DOS(u)*FermiDirac(u,sv.Tel,μ,gc)
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

function Cel_nonlinear(sv,gc,Tel,μ)::Float64
    int(u,p)=ForwardDiff.derivative(x->(exp((u-μ)/(gc.kB*x))+1)^-1,Tel)*sv.DOS(u)*u
    prob=IntegralProblem(int,sv.lb,sv.ub)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-10,abstol=1e-10)
    return abs(sol.u)::Float64
end

function Cel_linear(sv,mp) ::Float64
    return mp.γ*sv.Tel
end

function Cph(mp,sv,gc) ::Float64
    int(u,p)=(u^4*exp(u))/((exp(u)-1)^2)
    prob = IntegralProblem(int,0.0,mp.θ/sv.Tph)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-5,abstol=1e-5)
    return (9*mp.n*sv.kB)*(sv.Tph/mp.θ)^3 *sol.u
end 

function dTeldE(μ,E,Tel)
    v=1/(kB*Tel)
    u=-exp((E-μ)/(kB*Tel))
    w=(exp.((E-μ)/(kB*Tel))+1)^2
    f=v*u/w
    if isnan(f)
        f=0.0
    end
    return f
end

function gep_nonlinear(sv,gc,mp)
    int(u,p)=sv.DOS(u)^2*-ForwardDiff.derivative(x->(exp((u-sv.μ)/(gc.kB*sv.Tel))+1)^-1,u)
    prob=IntegralProblem(int,sv.lb,sv.ub)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-5,abstol=1e-5)
    sv.g= sol.u*(pi*gc.kB*mp.λ)/(sv.DOS(sv.μ)*gc.hbar)
end

function LaserStrength(tprime,l,mp)
    return mp.α*sqrt(4*log(2)/pi)*(l.ϕ/l.FWHM)*exp(-4*log(2)*((tprime-l.delay)/l.FWHM)^2)
end

function dTeldt(mp,sv,gc,l,sim,t)
    if sim.Cel=="L"
        Ce=Cel_linear(sv,mp)
    elseif Input["Cel"]=="NL"
        Ce=Cel_nonlinear(sv,gc,Tel,μ)
    end
    if sim.ElecPhon=="NL"
        sv.g=gep_nonlinear(sv,gc,mp)
    end
    if sim.vers==ETTM
        particle=particleconstant(DOS,μ,Tel)
        U=totalUee(sv,gc,l,mp,particle,τep,t)

    Teldt=(U)/Ce
    return Teldt
end

function dTphdt(Tel,Tph,t,g,μ,Input,DOS,τ)
    G=GTelTph(g,Tel,Tph,μ,Input,DOS)
    U=totalUep(g,Tph,τ,μ,Tel,t,DOS,Input)
    Tphdt=(G+U)/Cph(Tph)
    return Tphdt
end

function Trajectory()
    

export DensityOfStates,NumberOfElectrons,ChemicalPotential

end
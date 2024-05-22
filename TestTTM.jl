using ModelingToolkit,DifferentialEquations,Plots,Symbolics,Latexify,LaTeXStrings
using ModelingToolkit: t_nounits as t, D_nounits as D 

abstract type Laser end
@kwdef struct Gaussian <: Laser 
    FWHM::Real
    Offset::Real
    Fluence::Real
    hv::Real
    spatial::Array{Float64}
    R::Real
end

function (::Gaussian)(lp)
    temp = sqrt(4*log(2)/pi)/lp.FWHM*exp(-4*log(2)*((t-(2*lp.FWHM)-lp.Offset)/lp.FWHM)^2)
    return temp
end

function dTel(;name)
    @parameters g γ 
    @variables Tph(t) S(t) Tel(t) 
    
    eqs = D(Tel) ~ (-g*(Tel-Tph) + S)/(γ*Tel)

    ODESystem(eqs,t;name)

end

function dTph(;name)
    @parameters g Cph
    @variables Tel(t) Tph(t)
    
    eqs = D(Tph) ~ g*(Tel-Tph)/Cph

    ODESystem(eqs,t;name)

end

function Laser()
    @parameters  ϕ ϵ
    lp=Gaussian(FWHM=5e-14,Offset=2e-13,Fluence=10.0,hv=3.1,spatial=[0],R=0.0)
    return  lp(lp)*lp.Fluence/ϵ
end

function dTel_factory(;name)
    @named dTeldt = dTel()
    laser=Laser()
    connections=[dTeldt.S ~ laser]
    compose(ODESystem(connections,t;name),dTeldt)
end

function main()
    @named Tph_eq = dTph() 
    @named Tel_eq=dTel_factory()

    connections=[Tel_eq.dTeldt.Tph ~ Tph_eq.Tph
                Tph_eq.Tel ~ Tel_eq.dTeldt.Tel]
    
    connected = compose(ODESystem(connections,t,name=:connected
                ,defaults=Pair{Num,Any}[Tel_eq.dTeldt.g => Tph_eq.g])
                ,Tel_eq,Tph_eq)
    connected_simp=structural_simplify(connected)
    u0=[Tel_eq.dTeldt.Tel => 300.0,
        Tph_eq.Tph => 300,
        Tel_eq.dTeldt.S => 0.0]
    p=[Tph_eq.g => 2.3e16,
      Tel_eq.dTeldt.γ => 71,
      Tph_eq.Cph => 2.5e6,
      Tel_eq.ϵ => 1.27e-8]
    prob=ODEProblem(connected_simp,u0,(0.0,5e-13),p)
    sol=solve(prob,AutoVern7(Rodas5()),abstol=1e-10,reltol=1e-10)
    return sol
end


test=main()
plot(test)
using ModelingToolkit,DifferentialEquations,Symbolics,Dierckx,Integrals
using ForwardDiff,Cubature,Roots
using ModelingToolkit: t_nounits as t, D_nounits as D

#
#The next 250 lines defines the different individual systems
#
function athem_factory(laser::Num,egl::Int;name)
    @variables Tel(t) feq(t)[1:egl] n(t)
    @parameters kB egrid[1:egl] μ hv DOS::Spline1D u0 FE τ
    @named dfneq = athem_template(egl)
    feq = FermiDirac(egrid,μ,Tel,kB)
    connections = [dfneq.Source .~ athem_excitation(egrid,feq,dfneq.fneq,hv,DOS)*laser
                   dfneq.neqelel .~ -athem_electronelectronrelaxation(dfneq.fneq,feq,egrid,μ,kB,u0,n,DOS).*flt_relaxation(τ,FE,μ,egrid,kB,Tel)]

    connections = Symbolics.scalarize(reduce(vcat, Symbolics.scalarize.(connections)))

    compose(ODESystem(connections,t;name),dfneq)
end

function athem_template(egl::Int;name)
    @variables (fneq(t))[1:egl] (Source(t))[1:egl] (neqelel(t))[1:egl]

    eqs = [D(fneq) .~ Source .+ neqelel]

    eqs = Symbolics.scalarize(reduce(vcat, Symbolics.scalarize.(eqs)))

    ODESystem(eqs,t;name)
end

function athem_excitation(egrid,feq,fneq,hv,DOS)
    ftot = fneq.+feq
    ftotspl=get_interpolate(egrid,ftot)

    Δfe = fgr_electron_generation(egrid,DOS,ftotspl,hv)
    Δfh = fgr_hole_generation(egrid,DOS,ftotspl,hv)

    pc_sf = get_noparticles(Δfe,DOS,egrid) / get_noparticles(Δfh,DOS,egrid)
    Δfshape = (Δfe.*pc_sf).-Δfh
    inten = get_internalenergy_grid(Δfshape,DOS,egrid)
    return Δfshape./inten
end

function fgr_hole_generation(egrid,DOS::Spline1D,ftotspl::Spline1D,hv::Real)
    return DOS(egrid.+hv).*ftotspl(egrid).*(1 .-ftotspl(egrid.+hv))
end
@register_array_symbolic fgr_hole_generation(egrid,DOS,ftotspl,hv) begin
    size=(length(egrid),)
    eltype=eltype(egrid)
end
function fgr_electron_generation(egrid,DOS::Spline1D,ftotspl::Spline1D,hv::Real)
    return DOS(egrid.-hv).*ftotspl(egrid.-hv).*(1 .-ftotspl(egrid))
end
@register_array_symbolic fgr_electron_generation(egrid,DOS,ftotspl,hv) begin
    size=(length(egrid),)
    eltype=eltype(egrid)
end

function athem_electronelectronrelaxation(fneq,feq,egrid,μ,kB,u0,n,DOS)
    ftot = feq .+ fneq
    ftot_spl = get_interpolate(egrid,ftot)
    goal = get_internalenergyspl(μ,ftot_spl,DOS,u0)
    rel_T,rel_μ = find_relaxed_eedistribution(goal,kB,egrid,n,DOS,u0)
    return (fneq .+ frel .- feq)
end

function find_relaxed_eedistribution(goal,kB,egrid,n,DOS,u0)
    f(u) = goal - relaxed_internalenergy(u,n,DOS,kB,egrid,u0)
    rel_Tel = solve(ZeroProblem(f,1000.0),Order1();atol=1e-4,rtol=1e-4)
    rel_μ = find_chemicalpotential(n,rel_Tel,0.0,DOS,kB)
    return FermiDirac(egrid,rel_μ,rel_Tel,kB)
end
@register_array_symbolic find_relaxed_eedistribution(goal,kB,egrid,n,DOS,u0) begin
    size=(length(egrid),)
    eltype=eltype(egrid)
end

function relaxed_internalenergy(T,n,DOS,kB,egrid,u0)
    cp = find_chemicalpotential(n,T,0.0,DOS,kB)
    u = get_internalenergyspl(cp,get_interpolate(egrid,FermiDirac(egrid,cp,T,kB)),DOS,u0)
    return u
end

function flt_relaxation(τ,FE,μ,egrid,kB,Tel)
    return τ*(FE+μ)^2 ./((egrid.-(FE+μ)).^2 .+(pi*kB*Tel)^2)
end

function dFDdE(kB::Float64,Tel::Float64,μ::Float64,E::Float64)::Real
    Numer=-exp((E-μ)/(kB*Tel))
    Denom=kB*Tel*(exp((E-μ)/(kB*Tel))+1)^2
    return Numer/Denom
end

function dFDdT(kB::Float64,Tel::Float64,μ::Float64,E::Float64)::Real
    Numer=(E-μ)*exp((E-μ)/(kB*Tel))
    Denom=kB*Tel^2*(exp((E-μ)/(kB*Tel))+1)^2
    return Numer/Denom
end

function dFDdμ(kB::Float64,Tel::Float64,μ::Float64,E::Float64)::Real
    Numer=exp((E-μ)/(kB*Tel))
    Denom=kB*Tel*(exp((E-μ)/(kB*Tel))+1)^2
    return Numer/Denom
end

FermiDirac(E::Real,μ::Union{Real,ForwardDiff.Dual},Tel::Real,kB::Real)= 1/(exp((E-μ)/(kB*Tel))+1)

FermiDirac(egrid,μ,Tel,kB) = 1 ./(exp.((egrid.-μ)./(kB*Tel)).+1)

function particle_change(egl;name)
    @variables n(t) relax_dis(t)[1:egl]
    @parameters egrid[1:egl] μ n0 DOS::Spline1D
    
    eqs = [D(n) ~ relaxedelectron_particlechange(relax_dis,egrid,DOS,μ,n0)]

    ODESystem(eqs,t;name)
end

function relaxedelectron_particlechange(relax_dis,egrid,DOS::Spline1D,μ::Real,n0::Real)
    dis_spl = get_interpolate(egrid,relax_dis)
    return get_noparticlesspl(μ,Dis,DOS,n0)
end
@register_symbolic relaxedelectron_particlechange(relax_dis::AbstractVector,egrid::AbstractVector,DOS::Spline1D,μ::Num,n0::Num)

function t_electron_factory(egl::Real;name)
    @variables Tel(t) relax_dis(t)[1:egl]
    @parameters μ DOS::Spline1D egrid[1:egl] u0
    @named dTel = athem_elec_template()
    connections =[dTel.Δu ~ relaxedelectron_internalenergy(relax_dis,egrid,DOS,μ,u0)]
    compose(ODESystem(connections,t;name,defaults=[dTel.DOS=>DOS,dTel.μ=>μ]),dTel)

end

function athem_elec_template(;name)
    @variables Tel(t) Δu(t) Δn(t)
    @parameters kB μ DOS::Spline1D

    eqs = [D(Tel) ~ 1/(c_T(μ,Tel,DOS,kB)*p_μ(μ,Tel,DOS,kB)-p_T(μ,Tel,DOS,kB)*c_μ(μ,Tel,DOS,kB))*(p_μ(μ,Tel,DOS,kB)*Δu-c_μ(μ,Tel,DOS,kB)*Δn)]

    ODESystem(eqs,t;name)
end

function relaxedelectron_internalenergy(relax_dis::AbstractVector,egrid,DOS::Spline1D,μ::Real,u0::Real)
    dis_spl = get_interpolate(egrid,relax_dis)
    return get_internalenergyspl(μ,dis_spl,DOS,u0)
end
@register_symbolic relaxedelectron_internalenergy(relax_dis::AbstractVector,egrid::AbstractVector,DOS::SymbolicUtils.BasicSymbolic{Spline1D},μ::Num,u0::Num)::Real

get_interpolate(xvals::AbstractVector,yvals::AbstractVector) = Spline1D(xvals,yvals,bc="nearest")
@register_symbolic get_interpolate(xvals,yvals)::Spline1D

function update_chempotAthEM!(integ,u,p)
    integ.p[p.μ] = find_chemicalpotential(integ.u[u.n],integ.u[u.Tel],integ.p[p.μ],integ.p[p.DOS],integ.p[p.kB])
end

function find_chemicalpotential(no_part::Real,Tel::Real,μ::Real,DOS::Spline1D,kB::Real)
    f(u) = no_part - get_thermalparticles(u,Tel,DOS,kB)
    return solve(ZeroProblem(f,μ),Order1();atol=1e-3,rtol=1e-3)
end

function get_thermalparticles(μ,Tel::Real,DOS::Spline1D,kB::Real)
    int(u,p) = FermiDirac(u,μ,Tel,kB)*DOS(u)
    return solve(IntegralProblem(int,(μ-20,μ+20)),HCubatureJL(initdiv=50);abstol=1e-6,reltol=1e-6).u
end

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
@register_symbolic get_noparticles(Dis::AbstractVector,DOS,egrid::AbstractVector)

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

function get_internalenergyspl(μ::Real,Dis::Spline1D,DOS::Spline1D,u0::Real)
    int_neg(u,p) = (Dis(u).-1)*DOS(u)*u
    uroot_neg = solve(IntegralProblem(int_neg,(-Inf,0.0)),CubatureJLh();reltol=1e-3,abstol=1e-3).u
    int_pos(u,p) = Dis(u)*DOS(u)*u
    uroot_pos = solve(IntegralProblem(int_pos,(0.0,μ+10)),CubatureJLh();reltol=1e-3,abstol=1e-3).u
    return uroot_neg+uroot_pos+u0
end
@register_symbolic get_internalenergyspl(μ,Dis,DOS,u0)

function get_internalenergy_grid(Dis,DOS,egrid)
    integrand = Dis.*DOS(egrid).*egrid
    prob = SampledIntegralProblem(integrand,egrid)
    return solve(prob,SimpsonsRule()).u
end
@register_symbolic get_internalenergy_grid(Dis::AbstractVector,DOS,egrid::AbstractVector)
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

#
#This section assembles the systems
#

function system_builder()
    @parameters R FWHM ϕ ϵ
    laser = (0.9394372786996513(1 - R)*exp((-2.772588722239781(t^2)) / (FWHM^2))*ϕ) / (FWHM*ϵ)
    egl=605
    @named neq = athem_factory(laser,egl)
    @named elec_temp = t_electron_factory(egl)
    @named partchange = particle_change(egl)

    connections = [elec_temp.relax_dis ~ partchange.relax_dis
                   elec_temp.relax_dis ~ neq.dfneq.neqelel
                   elec_temp.dTel.Tel ~ neq.Tel
                   neq.n ~ partchange.n
                   elec_temp.dTel.Δn ~ D(partchange.n)]
    connections = Symbolics.scalarize(reduce(vcat, Symbolics.scalarize.(connections)))
    default_params = [partchange.egrid => elec_temp.egrid
                neq.egrid => elec_temp.egrid
                partchange.μ => elec_temp.μ
                neq.μ => elec_temp.μ
                neq.DOS => elec_temp.DOS
                partchange.DOS => elec_temp.DOS
                neq.kB => elec_temp.dTel.kB]
    events = (t>=-1e5) => (update_chempotAthEM!,[elec_temp.dTel.Tel=>:Tel,partchange.n => :n],
               [elec_temp.μ=>:μ,elec_temp.dTel.kB=>:kB,elec_temp.DOS=>:DOS],[elec_temp.μ])
    connected = compose(ODESystem(connections,t,name=:connected,defaults=default_params,discrete_events=events),neq,elec_temp,partchange)
    connected_sys = structural_simplify(connected)
    return connected_sys
end

system_builder()
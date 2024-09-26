
using ModelingToolkit,DifferentialEquations,Symbolics,Dierckx,Integrals,DelimitedFiles
using ForwardDiff,Cubature,Roots,StaticArrays,BenchmarkTools, JuliaSimCompiler
using Symbolics: scalarize
using ModelingToolkit: t_nounits as t, D_nounits as D

include("SimulationVariables.jl")
include("SimulationSetup.jl")
include("Lasers.jl")

#
#The next 250 lines defines the different individual systems
#
function athem_factory(laser::Num,egl::Int,DOS::Spline1D;name)
    @variables Tel(t) (feq(t))[1:egl] n(t) (fneq(t))[1:egl] (Source(t))[1:egl] (neqelel(t))[1:egl]
    @parameters μ hv egrid[1:egl] u0 n0 FE τ kB

    eqs = [scalarize(D.(fneq) .~ Source .+ neqelel)
           scalarize(feq .~ FermiDirac(egrid,μ,Tel,kB))
           scalarize(Source .~ athem_excitation(egrid,feq.+fneq,hv,DOS,μ,u0,n0).*laser)
           scalarize(neqelel .~ -athem_electronelectronrelaxation(fneq,feq,egrid,μ,kB,u0,n,DOS).*flt_relaxation(τ,FE,μ,egrid,kB,Tel))]

    ODESystem(eqs, t; name)
end

function athem_excitation(egrid,ftot,hv,DOS::Spline1D,μ,u0,n0)
    Δfe = fgr_electron_generation(egrid,DOS,ftot,hv)
    Δfh = fgr_hole_generation(egrid,DOS,ftot,hv)
    pc_sf = (get_noparticlesspl(μ,egrid,Δfe,DOS,n0) / get_noparticlesspl(μ,egrid,Δfh,DOS,n0))
    Δfshape = (Δfe*pc_sf).-Δfh
    inten = get_internalenergyspl(μ,egrid,Δfshape,DOS,u0)
    return Δfshape./inten
end
@register_array_symbolic athem_excitation(egrid::AbstractVector,ftot::AbstractVector,hv::Num,DOS::Spline1D,μ::Num,u0::Num,n0::Num) begin
    size = (length(ftot),)
    eltype=eltype(ftot)
end

function fgr_hole_generation(egrid,DOS::Spline1D,ftot,hv::Real)
    ftotspl=get_interpolate(egrid,ftot)
    return DOS(egrid.+hv).*ftotspl(egrid).*(1 .-ftotspl(egrid.+hv))
end
#= @register_array_symbolic fgr_hole_generation(egrid::AbstractVector,DOS::Spline1D,ftot::AbstractVector,hv::Num) begin
    size=(length(ftot),)
    eltype=eltype(ftot)
end =#
function fgr_electron_generation(egrid,DOS::Spline1D,ftot,hv::Real)
    ftotspl=get_interpolate(egrid,ftot)
    return DOS(egrid.-hv).*ftotspl(egrid.-hv).*(1 .-ftotspl(egrid))
end
#= @register_array_symbolic fgr_electron_generation(egrid::AbstractVector,DOS::Spline1D,ftot::AbstractVector,hv::Num) begin
    size=(length(ftot),)
    eltype=eltype(ftot)
end =#

function athem_electronelectronrelaxation(fneq,feq,egrid,μ,kB,u0,n,DOS::Spline1D)
    goal = get_internalenergyspl(μ,egrid,feq.+fneq,DOS,u0)
    frel = find_relaxed_eedistribution(goal,kB,egrid,n,DOS,u0)
    return fneq .+ frel .- feq
end
@register_array_symbolic athem_electronelectronrelaxation(fneq::AbstractVector,feq::AbstractVector,egrid::AbstractVector,μ::Num,kB::Num,u0::Num,n::Num,DOS::Spline1D) begin
    size = (length(fneq),)
    eltype = eltype(fneq)
end

function find_relaxed_eedistribution(goal,kB,egrid,n,DOS,u0)
    f(u) = goal - relaxed_internalenergy(u,n,DOS,kB,egrid,u0)
    rel_Tel = solve(ZeroProblem(f,1000.0),Order1();atol=1e-4,rtol=1e-4)
    rel_μ = find_chemicalpotential(n,rel_Tel,0.0,DOS,kB)
    return FermiDirac(egrid,rel_μ,rel_Tel,kB)
end
#= @register_array_symbolic find_relaxed_eedistribution(goal::Num,kB::Num,egrid::AbstractVector,n::Num,DOS::Spline1D,u0::Num) begin
    size=(length(egrid),)
    eltype=eltype(egrid)
end =#

function relaxed_internalenergy(T,n,DOS,kB,egrid,u0)
    cp = find_chemicalpotential(n,T,0.0,DOS,kB)
    u = get_internalenergyspl(μ,egrid,FermiDirac(egrid,μ,Tel,kB),DOS,u0)
    return u
end

function flt_relaxation(τ,FE,μ,egrid,kB,Tel)
    return τ*(FE+μ)^2 ./((egrid.-(FE+μ)).^2 .+(pi*kB*Tel)^2)
end

function dFDdT(kB::Float64,Tel::Float64,μ::Float64,E::Float64)
    Numer=(E-μ)*exp((E-μ)/(kB*Tel))
    Denom=kB*Tel^2*(exp((E-μ)/(kB*Tel))+1)^2
    return Numer/Denom
end

function dFDdμ(kB::Float64,Tel::Float64,μ::Float64,E::Float64)
    Numer=exp((E-μ)/(kB*Tel))
    Denom=kB*Tel*(exp((E-μ)/(kB*Tel))+1)^2
    return Numer/Denom
end

FermiDirac(egrid,μ,Tel,kB) = 1 ./(exp.((egrid.-μ)./(kB*Tel)).+1)

function particle_change(egl,DOS::Spline1D;name)
    @variables n(t) (relax_dis(t))[1:egl]
    @parameters μ egrid[1:egl] n0

    eqs = [D(n) ~ get_noparticlesspl(μ,egrid,relax_dis,DOS,n0)]

    ODESystem(eqs, t; name)
end

function t_electron_factory(egl::Real,DOS::Spline1D;name)
    @variables (relax_dis(t))[1:egl] Tel(t) Δu(t) Δn(t)
    @parameters μ egrid[1:egl] kB u0

    eqs = [D(Tel) ~ 1/(c_T(μ,Tel,DOS,kB)*p_μ(μ,Tel,DOS,kB)-p_T(μ,Tel,DOS,kB)*c_μ(μ,Tel,DOS,kB))*(p_μ(μ,Tel,DOS,kB)*Δu-c_μ(μ,Tel,DOS,kB)*Δn)
           Δu ~ get_internalenergyspl(μ,egrid,relax_dis,DOS,u0)]
    
    ODESystem(eqs,t; name)
end

get_interpolate(xvals,yvals) = Spline1D(xvals,yvals,bc="nearest")

function get_noparticlesspl(μ::Real,xvals,yvals,DOS::Spline1D,n0::Real)
    Dis=get_interpolate(xvals,yvals)
    int_neg(u,p) = (Dis(u).-1)*DOS(u)
    uroot_neg = solve(IntegralProblem(int_neg,(-Inf,0.0)),CubatureJLh();reltol=1e-3,abstol=1e-3).u
    int_pos(u,p) = Dis(u)*DOS(u)
    uroot_pos = solve(IntegralProblem(int_pos,(0.0,μ+10)),CubatureJLh();reltol=1e-3,abstol=1e-3).u
    return uroot_neg+uroot_pos+n0
end
@register_symbolic get_noparticlesspl(μ,xvals::AbstractVector,yvals::AbstractVector,DOS::Spline1D,n0::Num)

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

function get_internalenergyspl(μ::Real,xvals,yvals,DOS::Spline1D,u0::Real)
    Dis=get_interpolate(xvals,yvals)
    int_neg(u,p) = (Dis(u).-1)*DOS(u)*u
    uroot_neg = solve(IntegralProblem(int_neg,(-Inf,0.0)),CubatureJLh();reltol=1e-3,abstol=1e-3).u
    int_pos(u,p) = Dis(u)*DOS(u)*u
    uroot_pos = solve(IntegralProblem(int_pos,(0.0,μ+10)),CubatureJLh();reltol=1e-3,abstol=1e-3).u
    return uroot_neg+uroot_pos+u0
end
@register_symbolic get_internalenergyspl(μ::Num,xvals::AbstractVector,yvals::AbstractVector,DOS::Spline1D,u0::Num)

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

function get_n0(DOS,μ)
    int(u,p) = DOS(u)
    prob=IntegralProblem(int,(-Inf,μ))
    return solve(prob,HCubatureJL(initdiv=100),reltol=1e-5,abstol=1e-5).u
end

function get_u0(DOS,μ)
    int(u,p) = DOS(u)*u
    prob=IntegralProblem(int,(-Inf,μ))
    return solve(prob,HCubatureJL(initdiv=100),reltol=1e-5,abstol=1e-5).u
end

#
#This section assembles the systems
#

function setup()
    las=define_laser_system(:Gaussian,fwhm=25,fluence=111,photon_en=3.1)
    sim = define_simulation_settings(nlchempot=true,nlelecphon=true,nlelecheat=true,noneqelec=true
    ,elecphonint=false,elecelecint=true,electemp=true,phonontemp=false)
    mp = define_material_parameters(las,extcof=12.7,gamma=4.4315e-22,debye=165,noatoms=59,plasma=2.1357,thermalcond=320.0,
    elecperatom=1,eleceffmass=1.1,dos="DOS/Au_DOS.dat",secmomspecfun=23e-6,elecphon=1.44e-7,ballistic=0.0,cph=0.015)
    cons=Constants(8.617e-5,0.6582)
    dim = Homogenous()
    return sim,mp,las,dim,cons
end

function system_builder()
    sim,mp,las,dim,cons = setup()
    laser = laser_factory(las,mp,dim)
    u0=get_n0(mp.DOS,mp.μ)
    n0=get_u0(mp.DOS,mp.μ)

    egl=length(mp.egrid)

    @named neq = athem_factory(laser,egl,mp.DOS)
    println("fneq done")
    @named elec_temp = t_electron_factory(egl,mp.DOS)
    println("Tel done")
    @named partchange = particle_change(egl,mp.DOS)
    println("n done")

    connections = [elec_temp.relax_dis ~ neq.neqelel
                   partchange.relax_dis ~ neq.neqelel
                   elec_temp.Tel ~ neq.Tel
                   neq.n ~ partchange.n
                   elec_temp.Δn ~ D(partchange.n)]
    default_params = [partchange.μ => elec_temp.μ
                neq.μ => elec_temp.μ
                neq.egrid => elec_temp.egrid
                partchange.egrid => elec_temp.egrid
                neq.u0 => elec_temp.u0
                partchange.n0 => neq.n0
                neq.kB => elec_temp.kB]
    sys = ODESystem(connections, t, defaults=default_params; name=:sys)
    println("Connections done")
    @named model = compose(sys,neq,elec_temp,partchange)
    println("Begin Simplification")
    ir_model = IRSystem(model)
    connected_sys = structural_simplify(ir_model,fully_determined=false)
    return connected_sys
end

println("Compilation Done")
sys = system_builder()
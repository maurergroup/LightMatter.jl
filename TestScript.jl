using ModelingToolkit,OrdinaryDiffEq,Plots,Symbolics,Dierckx,DelimitedFiles,Integrals
using BenchmarkTools,ForwardDiff,StaticArrays,IfElse,Cubature,SimpleDiffEq,Roots
using ModelingToolkit: t_nounits as t, D_nounits as D

#= include("SymbolicsInterpolation.jl")
include("SimulationVariables.jl")
include("SimulationSetup.jl")
include("Lasers.jl")
include("ElectronDistribution.jl") =#

println("Compiled functions")

function setup()
    las=define_laser_system(:Gaussian,fwhm=25,fluence=111,photon_en=3.1)
    sim = define_simulation_settings(nlchempot=false,nlelecphon=true,nlelecheat=true,noneqelec=true                                                                                                          
    ,elecphonint=false,elecelecint=false,electemp=false,phonontemp=false)

    mp = define_material_parameters(las,extcof=12.7,gamma=4.4315e-22,debye=165,noatoms=59,plasma=2.1357,thermalcond=320.0,
    elecperatom=1,eleceffmass=1.1,dos="DOS/Au_DOS.dat",secmomspecfun=23e-6,elecphon=1.44e-7,ballistic=0.0,cph=0.015)
    cons=Constants(8.617e-5,0.6582)
    dim = Homogenous()
    return sim,mp,las,dim,cons
end


function eulermethod(time,egrid,source,hv,fneq,feq,DOS,dt)
    input = substitute(source,Dict([ϕ=>243,FWHM=>50.0,t=>time,ϵ=>12.7,R=>0.0]))
    Δ = athem_excitation(egrid,feq,fneq,hv,DOS)
    return Δ*dt*input.val
end

function euler_loop(trange,mp,laser,fneq,feq)
    for time in trange
        fneq .+= eulermethod(time,mp.egrid,laser,3.1,fneq,feq,mp.DOS,0.1)
    end
    return fneq
end

function main()
    sim,mp,las,dim,cons=setup()
    egl = length(mp.egrid)
    laser=laser_factory(las,dim)
    #= feq = 1 ./(exp.((mp.egrid.-mp.μ)./(cons.kB*300.0)).+1)
    fneq = zeros(egl)
    trange = range(-250.0,250.0,step=0.1)

    #= fneq = euler_loop(trange,mp,laser,fneq,feq)
    return fneq =#
    @btime fneq = euler_loop($trange,$mp,$laser,$fneq,$feq) =#
    
    @named neq = athem_factory(sim,mp.DOS,laser,egl)
    return neq
    #= simp = structural_simplify(neq)
    u0=[simp.dfneq.fneq=> zeros(egl)]
    p=[simp.μ=>mp.μ,
    simp.egrid=>mp.egrid,
    simp.FWHM=>las.FWHM,
    simp.kB=>cons.kB,
    simp.hv=>las.hv,
    simp.ϵ=>mp.ϵ,
    simp.ϕ=>las.ϕ,
    simp.Tel=>300.0,
    simp.R=>las.R]
    prob=ODEProblem(simp,u0,(-100.0,100.0),p)
    println("Start Solving")
    sol =solve(prob,Tsit5();dtmin=0.1,abstol=1e-7,reltol=1e-7,dense=false)
    return sol =#
end

sim,mp,las,dim,cons = setup()
laser=laser_factory(las,dim)
egl = length(mp.egrid)

@named test = athem_factory(sim,mp.DOS,laser,egl)
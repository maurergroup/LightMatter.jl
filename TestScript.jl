using ModelingToolkit,OrdinaryDiffEq,Plots,Symbolics,Dierckx,DelimitedFiles,Integrals
using BenchmarkTools,ForwardDiff,StaticArrays,IfElse,Cubature,SimpleDiffEq,Roots
using ModelingToolkit: t_nounits as t, D_nounits as D

include("SymbolicsInterpolation.jl")
include("SimulationVariables.jl")
include("SimulationSetup.jl")
include("Lasers.jl")
include("ElectronTemperature.jl")
include("PhononTemperature.jl")
include("ElectronDistribution.jl")
include("SystemBuilder.jl")

println("Compiled functions")

function setup()
    las=define_laser_system(:Gaussian,fwhm=50,fluence=243,photon_en=3.1)
    sim = define_simulation_settings(nlchempot=false,nlelecphon=true,nlelecheat=true,noneqelec=true                                                                                                          
    ,elecphonint=false,elecelecint=false,electemp=false,phonontemp=false)
    mp = define_material_parameters(las,extcof=12.7,gamma=4.4315e-22,debye=165,noatoms=59,plasma=2.1357,thermalcond=320.0,
    elecperatom=1,eleceffmass=1.1,dos="DOS/Au_DOS.dat",secmomspecfun=23e-6,elecphon=1.44e-7,ballistic=0.0,cph=0.015)
    cons=Constants(8.617e-5,0.6582)
    dim = Homogenous()
    return sim,mp,las,dim,cons
end

function main()
    sim,mp,las,dim,cons=setup()
    egl = length(mp.egrid)
    laser=laser_factory(las,dim)

    @named neq = athem_factory(sim,mp.DOS,laser,egl)
    simp = structural_simplify(neq)
    u0=[dfneq.fneq=> @SVector zeros(egl)]
    p=[μ=>mp.μ,
    egrid=>mp.egrid,
    FWHM=>las.FWHM,
    kB=>cons.kB,
    hv=>las.hv,
    ϵ=>mp.ϵ,
    ϕ=>las.ϕ,
    Tel=>300.0,
    R=>las.R]
    prob=ODEProblem(simp,u0,(-250.0,250.0),p)
    println("Start Solving")
    sol=solve(prob,SimpleATsit5();dt=0.1,abstol=1e-3,reltol=1e-3)
    return sol,simp
end
sol,simp=main()
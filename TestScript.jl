using ModelingToolkit,DifferentialEquations,Plots,Symbolics,Dierckx,DelimitedFiles,Integrals
using BenchmarkTools,ForwardDiff,StaticArrays,IfElse,Cubature
using ModelingToolkit: t_nounits as t, D_nounits as D
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

include("SymbolicsInterpolation.jl")
include("SimulationVariables.jl")
include("SimulationSetup.jl")
include("Lasers.jl")
include("ElectronTemperature.jl")
include("PhononTemperature.jl")
include("ElectronDistribution.jl")
include("SystemBuilder.jl")

function setup()
    las=define_laser_system(:Gaussian,fwhm=50,fluence=243,photon_en=3.1)
    sim = define_simulation_settings(nlchempot=true,nlelecphon=true,nlelecheat=true,noneqelec=false,elecphonint=true,elecelecint=false)
    mp = define_material_parameters(las,extcof=12.7,gamma=4.4315e-22,debye=165,noatoms=59,plasma=2.1357,thermalcond=320.0,
    elecperatom=1,eleceffmass=1.1,dos="DOS/Au_DOS.dat",secmomspecfun=23e-6,elecphon=1.44e-7,ballistic=0.0,cph=0.015)
    cons=Constants(8.617e-5,0.6582)
    dim = Homogenous()
    return sim,mp,las,dim,cons
end

function main()
    sim,mp,las,dim,cons=setup()
    egl = length(mp.egrid)
    fneq=zeros(egl)
    laser=laser_factory(las,dim)

    @named test_eq = athem_factory(mp.DOS,laser,egl)
    simp = structural_simplify(test_eq)

    u0=[simp.dfneq.fneq=>@SVector zeros(egl)]
    p=[simp.μ=>mp.μ,
    simp.egrid=>mp.egrid,
    simp.FWHM=>las.FWHM,
    simp.kB=>cons.kB,
    simp.hv=>las.hv,
    simp.ϵ=>mp.ϵ,
    simp.ϕ=>las.ϕ,
    simp.Tel=>300.0,
    simp.R=>las.R]
    prob=ODEProblem(simp,u0,(-250.0,250.0),p)
    @btime solve(prob,RK4();abstol=1e-5,reltol=1e-5,dtmin=0.1,save_everystep=false)
end

main()
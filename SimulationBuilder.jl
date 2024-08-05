using ModelingToolkit,DifferentialEquations,Plots,Symbolics,Dierckx,DelimitedFiles,Integrals
using Unitful,BenchmarkTools,ForwardDiff,StaticArrays,IfElse
using ModelingToolkit: t_nounits as t, D_nounits as D 

include("SymbolicsInterpolation.jl")
include("SimulationVariables.jl")
include("SimulationSetup.jl")
include("Lasers.jl")
include("ElectronTemperature.jl")
include("PhononTemperature.jl")
include("ElectronDistribution.jl")


function setup()
    las=define_laser_system(:Gaussian,fwhm=50,fluence=243,photon_en=3.1)
    sim = define_simulation_settings(nlchempot=true,nlelecphon=true,nlelecheat=true,noneqelec=false,elecphonint=true,elecelecint=false)
    mp = define_material_parameters(las,extcof=12.7,gamma=4.4315e-22,debye=165,noatoms=59,plasma=2.1357,thermalcond=320.0,
    elecperatom=1,eleceffmass=1.1,dos="DOS/Au_DOS.dat",secmomspecfun=23e-6,elecphon=1.44e-7,ballistic=0.0,cph=0.015)
    cons=Constants(8.617e-5,0.6582)
    dim = Homogenous()
    laser=laser_factory(las,dim)
    return sim,mp,las,laser,dim,cons
end

function run_dynamics(connected_eq,u0,tspan)

    prob=ODEProblem(connected_eq,u0,tspan,p)
    sol=solve(prob,Tsit5();abstol=1e-3,reltol=1e-3)
    return sol
end

function main()
    sim,mp,las,laser,dim=setup()
    connected_eq,Tel_eq,Tph_eq=equation_builder(sim,mp,laser)
    sol=run_dynamics(connected_eq,Tel_eq,Tph_eq,las,mp)
    p1=plot(sol,title="Temperature Profile",ylabel="Temperature / K",label=["Tel" "Tph"]
    ,linewidth=3,titlefont=12)
    p2=plot(sol,idxs=[Tel_eq.dTel.Source],title="Laser Profile",ylabel="Power",
    xlabel="Time / fs",label="Rectangular",linewidth=3,titlefont=12)
    plot(p1,p2,layout=[2,1],plot_title="Rectangular Laser Results")
end

#main()
sim,mp,las,laser,dim,cons=setup()
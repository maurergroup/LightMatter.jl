using DataInterpolations,DelimitedFiles,Integrals,Plots,Roots,RecursiveArrayTools,OrdinaryDiffEq,JLD2,FastGaussQuadrature,HDF5
include("SimulationSetup.jl")
include("Lasers.jl")
include("ElectronicTemperature.jl")
include("PhononicTemperature.jl")
include("SimulationVariables.jl")
include("ElectronicDistribution.jl")
include("SimulationConfigurations.jl")
include("SystemBuilder.jl")

function setup()
    las=define_laser_system(:Gaussian,fwhm=150,fluence=124.8,photon_en=2.18)
    sim = define_simulation_settings(nlelecphon=true,nlelecheat=true,noneqelec=true
    ,elecphonint=true,elecelecint=true,electemp=true,phonontemp=true)
    mp = define_material_parameters(las,extcof=14.9,gamma=6.117e-7,debye=343,noatoms=85,plasma=13.4,thermalcond=0.0025,
    elecperatom=1,eleceffmass=1.01,dos="DOS/Cu_DOS.dat",secmomspecfun=29e-6,elecphon=6.24e-7,ballistic=0.0,cph=0.015,τf=24.9)
    cons=Constants(8.617e-5,0.6582)
    dim = define_sim_dimensions(Dimension=1,Lengths=400,spacing=2)#define_sim_dimensions(Dimension=0)#
    return sim,mp,las,dim,cons
end

key_list = function_builder()
initialtemps=Dict("Tel"=>300.0,"Tph"=>300.0)
tspan=(-450, 1000.0)
#= sol = run_simulation(key_list,initialtemps,tspan)
@save "Test.jld2" sol =#
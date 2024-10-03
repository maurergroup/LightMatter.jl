using DataInterpolations,DelimitedFiles,Integrals,Plots,Roots,RecursiveArrayTools,OrdinaryDiffEq,DataFrames
include("SimulationSetup.jl")
include("Lasers.jl")
include("ElectronicTemperature.jl")
include("PhononicTemperature.jl")
include("SimulationVariables.jl")
include("ElectronicDistribution.jl")
include("SimulationConfigurations.jl")
include("SystemBuilder.jl")

function setup()
    las=define_laser_system(:Gaussian,fwhm=150,fluence=124.8,photon_en=1.55)
    sim = define_simulation_settings(nlchempot=true,nlelecphon=true,nlelecheat=true,noneqelec=false
    ,elecphonint=true,elecelecint=false,electemp=true,phonontemp=true)
    mp = define_material_parameters(las,extcof=14.9,gamma=4.4315e-22,debye=343,noatoms=85,plasma=13.4,thermalcond=0.0025,
    elecperatom=1,eleceffmass=1.01,dos="DOS/Cu_DOS.dat",secmomspecfun=29e-6,elecphon=6.24e-7,ballistic=0.0,cph=0.015,τf=18.55)
    cons=Constants(8.617e-5,0.6582)
    dim = define_sim_dimensions(Dimension=0)#define_sim_dimensions(Dimension=1,Lengths=400,spacing=2)#
    return sim,mp,las,dim,cons
end

key_list = function_builder()
initialtemps=Dict("Tel"=>300.0,"Tph"=>300.0)
tspan=(-450.0,550.0)
sol = run_simulation(key_list,initialtemps,tspan)
df - DataFrame(sol)
CSV.write("Test.csv", df)
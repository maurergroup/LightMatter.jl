using StaticArrays,Dierckx,ForwardDiff,DelimitedFiles,Integrals,OrdinaryDiffEq,Plots

include("MetaSimulationSetup.jl")
include("MetaLasers.jl")
include("MetaElectronicTemperature.jl")
include("MetaPhononicTemperature.jl")
include("MetaSimulationVariables.jl")


function setup()
    sim = define_simulation_settings(nlchempot=true,nlelecphon=true,nlelecheat=true,noneqelec=false,elecphonint=true)
    mp = define_material_parameters(extcof=12.7,gamma=4.4315e-22,debye=165,noatoms=59,plasma=2.1357,thermalcond=320.0,
    elecperatom=1,eleceffmass=1.1,dos="DOS/Au_DOS.dat",secmomspecfun=23e-6,elecphon=1.44e-7,ballistic=0.0,cph=0.015)
    las=define_laser_system(:Gaussian,fwhm=50,fluence=62.42,photon_en=3.1)
    laser=laser_factory(las,mp)
    return sim,mp,las,laser
end

function simulation(du,u,p,t)
    du[1] = Tel_eq(u[1],u[2],0.0,p[1],p[2],p[3],t)
    du[2] = Tph_eq(u[1],u[2],0.0,p[2],p[3],t)
end

function run_dynamics(las,mp,cons)
    u0 = [300.0,300.0]
    p=(las,mp,cons)
    prob=ODEProblem(simulation,u0,(-200,200),p)
    sol = solve(prob,Tsit5();abstol=1e-3,reltol=1e-3)
    return sol
end


sim,mp,las,laser = setup()
cons=Constants(8.617e-5,0.6582)
dim = Homogenous()
Tel_expr = electrontemperature_factory(mp,cons,sim,laser,dim)
Tph_expr = phonontemperature_factory(mp,cons,sim)
Tel_eq = eval(:((Tel,Tph,μ,las,mp,cons,t) -> $Tel_expr))
Tph_eq = eval(:((Tel,Tph,μ,mp,cons,t) -> $Tph_expr))
@time sol = run_dynamics(las,mp,cons)


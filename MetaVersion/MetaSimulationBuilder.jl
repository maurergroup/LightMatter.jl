using StaticArrays,Dierckx,ForwardDiff,DelimitedFiles,Integrals,OrdinaryDiffEq,Plots,NonlinearSolve

include("MetaSimulationSetup.jl")
include("MetaLasers.jl")
include("MetaElectronicTemperature.jl")
include("MetaPhononicTemperature.jl")
include("MetaSimulationVariables.jl")
include("MetaElectronicDistribution.jl")


function setup()
    las=define_laser_system(:Gaussian,fwhm=50,fluence=62.42,photon_en=3.1)
    sim = define_simulation_settings(nlchempot=true,nlelecphon=true,nlelecheat=true,noneqelec=true,elecphonint=false,elecelecint=true)
    mp = define_material_parameters(las,extcof=12.7,gamma=4.4315e-22,debye=165,noatoms=59,plasma=2.1357,thermalcond=320.0,
    elecperatom=1,eleceffmass=1.1,dos="DOS/Au_DOS.dat",secmomspecfun=23e-6,elecphon=1.44e-7,ballistic=0.0,cph=0.015)
    cons=Constants(8.617e-5,0.6582)
    dim = Homogenous()
    laser=laser_factory(las,dim)
    return sim,mp,las,laser,cons,dim
end

#= function simulation(du,u,p,t)
    du[1] = Tel_eq(u[1],u[2],0.0,p[1],p[2],p[3],t)
    du[2] = Tph_eq(u[1],u[2],0.0,p[2],p[3],t)
end =#

function run_dynamics(p)
    u0 = zeros(length(mp.egrid))
    simulation(u,p,t) = neq_eq(u,300.0,300.0,0.0,mp,cons,las,t,noe)
    prob=ODEProblem(simulation,u0,(-200,200),p)
    sol = solve(prob,Tsit5();abstol=1e-5,reltol=1e-5,saveat=1.0)
    return sol
end


sim,mp,las,laser,cons,dim = setup()
#= Tel_expr = electrontemperature_factory(sim,laser,dim)
Tph_expr = phonontemperature_factory(sim)
Tel_eq = eval(:((Tel,Tph,μ,las,mp,cons,t) -> $Tel_expr))
Tph_eq = eval(:((Tel,Tph,μ,mp,cons,t) -> $Tph_expr)) =#
noe = get_thermalparticles(0.0,1e-16,mp.DOS,cons.kB)
neq_expr = athemdistribution_factory(sim,laser)
neq_eq = eval(:((fneq,Tel,Tph,μ,mp,cons,las,t,no_part) -> $neq_expr))
p=(las,mp,cons,noe)
sol = run_dynamics(p)


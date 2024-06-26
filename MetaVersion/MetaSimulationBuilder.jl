using StaticArrays,Dierckx,ForwardDiff,DelimitedFiles,Integrals,OrdinaryDiffEq,Plots,Roots,Cubature,RecursiveArrayTools,NonlinearSolve

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

function simulation(du,u,p,t)
    println(t)
    du.x[1][1] = Tel_eq(u.x[1][1],u.x[2][1],t,p[1],p[2],p[3],p[4],u.x[3])
    du.x[2][1] = Tph_eq(u.x[1][1],u.x[2][1],t,p[1],p[3],p[4])
    du.x[3][:] = neq_eq(u.x[3],u.x[1][1],u.x[2][1],t,p[1],p[5],p[2],p[3],p[4])
end

function condition(u,t,integrator)
    true
end

function update_chempot!(integrator)
    integrator.p[1] = find_chemicalpotential(integrator.p[5],integrator.u.x[1][1],integrator.p[1],integrator.p[3].DOS,integrator.p[4].kB)
end

function run_dynamics(p)
    u0 = ArrayPartition([300.0],[300.0],zeros(length(mp.egrid)))
    prob=ODEProblem(simulation,u0,(-200,200),p)
    if sim.ParameterApprox.ChemicalPotential == true
        chempot = DiscreteCallback(condition,update_chempot!)
        sol = solve(prob,Tsit5();callback=chempot,abstol=1e-3,reltol=1e-3,saveat=1.0)
    elseif sim.ParameterApprox.ChemicalPotential == false
        sol = solve(prob,Tsit5();abstol=1e-3,reltol=1e-3,saveat=1.0)
    end
    return sol
end


sim,mp,las,laser,cons,dim = setup()
Tel_expr = electrontemperature_factory(sim,laser,dim)
Tph_expr = phonontemperature_factory(sim)
Tel_eq = eval(:((Tel,Tph,t,μ,las,mp,cons,fneq) -> $Tel_expr))
Tph_eq = eval(:((Tel,Tph,t,μ,mp,cons) -> $Tph_expr))
noe = get_thermalparticles(0.0,1e-16,mp.DOS,cons.kB)
neq_expr = athemdistribution_factory(sim,laser)
neq_eq = eval(:((fneq,Tel,Tph,t,μ,no_part,las,mp,cons) -> $neq_expr))
p=[0.0,las,mp,cons,noe]
sol = run_dynamics(p)


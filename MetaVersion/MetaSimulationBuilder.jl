using StaticArrays,Dierckx,ForwardDiff,DelimitedFiles,Integrals,OrdinaryDiffEq,Plots,Roots,Cubature,RecursiveArrayTools
include("MetaSimulationSetup.jl")
include("MetaLasers.jl")
include("MetaElectronicTemperature.jl")
include("MetaPhononicTemperature.jl")
include("MetaSimulationVariables.jl")
include("MetaElectronicDistribution.jl")

function setup()
    las=define_laser_system(:Gaussian,fwhm=50,fluence=243,photon_en=3.1)
    sim = define_simulation_settings(nlchempot=true,nlelecphon=true,nlelecheat=true,noneqelec=true,elecphonint=false,elecelecint=false)
    mp = define_material_parameters(las,extcof=12.7,gamma=4.4315e-22,debye=165,noatoms=59,plasma=2.1357,thermalcond=320.0,
    elecperatom=1,eleceffmass=1.1,dos="DOS/Au_DOS.dat",secmomspecfun=23e-6,elecphon=1.44e-7,ballistic=0.0,cph=0.015)
    cons=Constants(8.617e-5,0.6582)
    dim = Homogenous()
    laser=laser_factory(las,dim)
    return sim,mp,las,laser,cons,dim
end

function generate_expressions(sim,laser,dim)
    exprs = Dict{String,Expr}()
    if sim.Systems.ElectronTemperature == true
        merge!(exprs,Dict("Tel" => electrontemperature_factory(sim,laser,dim)))
    end
    if sim.Systems.PhononTemperature == true
        merge!(exprs,Dict("Tph" => phonontemperature_factory(sim)))
    end
    if sim.Systems.NonEqElectrons == true
        merge!(exprs,Dict("fneq" => athemdistribution_factory(sim,laser)))
        merge!(exprs,Dict("noe" => neqelectron_electronparticlechange()))
    end
    return exprs
end

function eval_expr(expr,vars)
    func = eval(Expr(:->, Expr(:parameters,keys(vars)...),expr))
    return Base.invokelatest(func;vars...)
end

function simulation(du,u,p,t)
    println(t)
    μ = find_chemicalpotential(u.x[4][1],u.x[1][1],p[3].FE,p[3].DOS,p[4].kB)
    vars = (Tel=u.x[1][1],Tph=u.x[2][1],fneq=u.x[3][:],t=t,noe=u.x[4],μ=μ,las=p[2],cons=p[4],mp=p[3])
    for i in 1:length(u.x)
        du.x[i][:] = eval_expr(p[1][i],vars)
    end
end

function run_dynamics(sim,p,u0)
    prob=ODEProblem(simulation,u0,(-150,200.0),p)
     sol = solve(prob,Vern7();callback=cbs,abstol=1e-5,reltol=1e-5,saveat=1.0,dtmin=0.05,force_dtmin=true)
    return sol
end

#= function main()
    #Build simulation with settings
    sim,mp,las,laser,cons,dim = setup()
    exprs = generate_expressions(sim,laser,dim)
    
    #= p=[exprs,las,mp,cons]
    u0 = ArrayPartition([300.0],[300.0],zeros(length(mp.egrid)),[get_thermalparticles(mp.FE,1e-16,mp.DOS,cons.kB)])

    sol=run_dynamics(sim,p,u0)
    return sol =#
end =#


#sol=main()

fneqout = sol(200.0,idxs=3:length(mp.egrid)+2)

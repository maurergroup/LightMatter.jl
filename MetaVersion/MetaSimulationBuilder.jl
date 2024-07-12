using StaticArrays,Dierckx,ForwardDiff,DelimitedFiles,Integrals,OrdinaryDiffEq,Plots,Roots,Cubature,RecursiveArrayTools,NonlinearSolve
include("MetaSimulationSetup.jl")
include("MetaLasers.jl")
include("MetaElectronicTemperature.jl")
include("MetaPhononicTemperature.jl")
include("MetaSimulationVariables.jl")
include("MetaElectronicDistribution.jl")

function eval_expr(expr,vars)
    func = eval(Expr(:->, Expr(:parameters,keys(vars)...),expr))
    return Base.invokelatest(func;vars...)
end

function setup()
    las=define_laser_system(:Gaussian,fwhm=50,fluence=243,photon_en=3.1)
    sim = define_simulation_settings(nlchempot=true,nlelecphon=true,nlelecheat=true,noneqelec=true,elecphonint=false,elecelecint=true)
    mp = define_material_parameters(las,extcof=12.7,gamma=4.4315e-22,debye=165,noatoms=59,plasma=2.1357,thermalcond=320.0,
    elecperatom=1,eleceffmass=1.1,dos="DOS/Au_DOS.dat",secmomspecfun=23e-6,elecphon=1.44e-7,ballistic=0.0,cph=0.015)
    cons=Constants(8.617e-5,0.6582)
    dim = Homogenous()
    laser=laser_factory(las,dim)
    return sim,mp,las,laser,cons,dim
end

#= function simulation(du,u,p,t)
    println(t)
    vars = (Tel=u.x[1][1],Tph=u.x[2][1],fneq=u.x[3][:],t=t,noe=p[5],μ=p[1],las=p[2],cons=p[4],mp=p[3])
    du.x[1][1] = eval_expr(p[6][1],vars)
    du.x[2][1] = eval_expr(p[6][2],vars)
    du.x[3][:] = eval_expr(p[6][3],vars)
end =#

function simulation(du,u,p,t)
    println(t)
    du.x[1][1] = Tel_eq(u.x[1][1],u.x[2][1],u.x[3][:],t,p[5],p[1],p[2],p[4],p[3])
    du.x[2][1] = Tph_eq(u.x[1][1],u.x[2][1],u.x[3][:],t,p[5],p[1],p[2],p[4],p[3])
    du.x[3][:] = neq_eq(u.x[1][1],u.x[2][1],u.x[3][:],t,p[5],p[1],p[2],p[4],p[3])
end

function condition(u,t,integrator)
    true
end

function update_chempot!(integrator)
    integrator.p[1] = find_chemicalpotential(integrator.p[5],integrator.u.x[1][1],integrator.p[1],integrator.p[3].DOS,integrator.p[4].kB)
    push!(integrator.p[6],integrator.p[1])
end

function update_noe!(integrator)
    vars = (Tel=integrator.u.x[1][1],cons=integrator.p[4],mp=integrator.p[3],μ=integrator.p[1],fneq=integrator.u.x[3][:])
    integrator.p[5] += eval_expr(neqelectron_electronparticlechange(),vars)
    push!(integrator.p[6][2],integrator.p[5])
end

function run_dynamics(sim,p)
    u0 = ArrayPartition([300.0],[300.0],zeros(length(p[3].egrid)))
    prob=ODEProblem(simulation,u0,(-150,500.0),p)
    if sim.ParameterApprox.ChemicalPotential == true
        chempot = DiscreteCallback(condition,update_chempot!)
        #nothermalelectron = DiscreteCallback(condition,update_noe!)
        #cbs = CallbackSet(chempot,nothermalelectron)
        sol = solve(prob,RK4();callback=chempot,abstol=1e-4,reltol=1e-4,saveat=1.0,dtmin=0.1,force_dtmin=true)
    elseif sim.ParameterApprox.ChemicalPotential == false
        sol = solve(prob,Tsit5();abstol=1e-3,reltol=1e-3,saveat=1.0)
    end
    return sol
end

#= function main()
    sim,mp,las,laser,cons,dim = setup()
    Tel_expr = electrontemperature_factory(sim,laser,dim)
    Tph_expr = phonontemperature_factory(sim)
    neq_expr = athemdistribution_factory(sim,laser)

    noe = get_thermalparticles(mp.FE,1e-16,mp.DOS,cons.kB)
    μ = find_chemicalpotential(noe,300.0,mp.FE,mp.DOS,cons.kB)
    exprs=[Tel_expr,Tph_expr,neq_expr]
    storage=[zeros(0),zeros(0)]
    p=[μ,las,mp,cons,noe,exprs,storage]
    sol = run_dynamics(sim,p)
    return sol
end

sol = main() =#

function RungeKutta(Tel,Tph,fneq,t,dt,noe,μ,las,cons,mp)
    Tel1 = dt*Tel_eq(Tel,Tph,fneq,t,noe,μ,las,cons,mp)
    Tph1 = dt*Tph_eq(Tel,Tph,fneq,t,noe,μ,las,cons,mp)
    fneq1 = dt*neq_eq(Tel,Tph,fneq,t,noe,μ,las,cons,mp)
    
    Tel2 = dt*Tel_eq(Tel+Tel1/2,Tph+Tph1/2,fneq+fneq1/2,t+dt/2,noe,μ,las,cons,mp)
    Tph2 = dt*Tph_eq(Tel+Tel1/2,Tph+Tph1/2,fneq+fneq1/2,t+dt/2,noe,μ,las,cons,mp)
    fneq2 = dt*neq_eq(Tel+Tel1/2,Tph+Tph1/2,fneq+fneq1/2,t+dt/2,noe,μ,las,cons,mp)

    Tel3 = dt*Tel_eq(Tel+Tel2/2,Tph+Tph2/2,fneq+fneq2/2,t+dt/2,noe,μ,las,cons,mp)
    Tph3 = dt*Tph_eq(Tel+Tel2/2,Tph+Tph2/2,fneq+fneq2/2,t+dt/2,noe,μ,las,cons,mp)
    fneq3 = dt*neq_eq(Tel+Tel2/2,Tph+Tph2/2,fneq+fneq2/2,t+dt/2,noe,μ,las,cons,mp)

    Tel4 = dt*Tel_eq(Tel+Tel3,Tph+Tph3,fneq+fneq3,t+dt,noe,μ,las,cons,mp)
    Tph4 = dt*Tph_eq(Tel+Tel3,Tph+Tph3,fneq+fneq3,t+dt,noe,μ,las,cons,mp)
    fneq4 = dt*neq_eq(Tel+Tel3,Tph+Tph3,fneq+fneq3,t+dt,noe,μ,las,cons,mp)

    Tel=Tel+(1/6*(Tel1+(2*Tel2)+(2*Tel3)+Tel4))
    Tph=Tph+(1/6*(Tph1+(2*Tph2)+(2*Tph3)+Tph4))
    fneq=fneq+(1/6*(fneq1+(2*fneq2)+(2*fneq3)+fneq4))

    μ = find_chemicalpotential(noe,Tel,μ,mp.DOS,cons.kB)
    vars = (Tel=Tel,cons=cons,mp=mp,μ=μ,fneq=fneq)
    noe += eval_expr(neqelectron_electronparticlechange(),vars)

    return Tel,Tph,fneq,μ,noe
end

sim,mp,las,laser,cons,dim = setup()
Tel_expr = electrontemperature_factory(sim,laser,dim)
Tph_expr = phonontemperature_factory(sim)
neq_expr = athemdistribution_factory(sim,laser)

Tel_eq = eval(:((Tel,Tph,fneq,t,noe,μ,las,cons,mp) -> $Tel_expr))
Tph_eq  = eval(:((Tel,Tph,fneq,t,noe,μ,las,cons,mp) -> $Tph_expr))
neq_eq = eval(:((Tel,Tph,fneq,t,noe,μ,las,cons,mp) -> $neq_expr))

Tel = 300.0
Tph= 300.0
fneq = zeros(length(mp.egrid))
noe = get_thermalparticles(mp.FE,1e-16,mp.DOS,cons.kB)
μ = find_chemicalpotential(noe,300.0,mp.FE,mp.DOS,cons.kB)
storage = []
p=[μ,las,mp,cons,noe,storage]
sol=run_dynamics(sim,p)

#= trange = range(-150.0,200.0,step=0.1)
Telout = zeros(length(trange))
fneqout = zeros(length(trange),length(mp.egrid))

for (i,t) in enumerate(trange)
    global Tel,Tph,fneq,noe,μ,las,cons,mp
    println(t)
    Tel += Tel_eq(Tel,Tph,fneq,t,noe,μ,las,cons,mp)
    fneq = fneq .+ neq_eq(Tel,Tph,fneq,t,noe,μ,las,cons,mp)
    μ = find_chemicalpotential(noe,Tel,mp.FE,mp.DOS,cons.kB)
    vars = (Tel=Tel,cons=cons,mp=mp,μ=μ,fneq=fneq)
    noe += eval_expr(neqelectron_electronparticlechange(),vars)

    Telout[i] = Tel
    fneqout[i,:]=fneq
end =#
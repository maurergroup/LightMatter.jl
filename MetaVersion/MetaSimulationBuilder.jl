using StaticArrays,Interpolations,ForwardDiff,DelimitedFiles,Integrals,OrdinaryDiffEq,Plots,Roots,Cubature,RecursiveArrayTools
include("MetaSimulationSetup.jl")
include("MetaLasers.jl")
include("MetaElectronicTemperature.jl")
include("MetaPhononicTemperature.jl")
include("MetaSimulationVariables.jl")
include("MetaElectronicDistribution.jl")

function setup()
    las=define_laser_system(:Gaussian,fwhm=25,fluence=111,photon_en=3.1)
    sim = define_simulation_settings(nlchempot=true,nlelecphon=true,nlelecheat=true,noneqelec=true
    ,elecphonint=false,elecelecint=true,electemp=true,phonontemp=false)
    mp = define_material_parameters(las,extcof=12.7,gamma=4.4315e-22,debye=165,noatoms=59,plasma=13.719,thermalcond=320.0,
    elecperatom=1,eleceffmass=1.1,dos="DOS/Au_DOS.dat",secmomspecfun=23e-6,elecphon=1.44e-7,ballistic=0.0,cph=0.015)
    cons=Constants(8.617e-5,0.6582)
    dim = Homogenous()
    return sim,mp,las,dim,cons
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
        if sim.Interactions.ElectronElectron == true
            merge!(exprs,Dict("noe" => athem_electronparticlechange()))
        end
    end
    return exprs,collect(keys(exprs))
end

function eval_expr(expr,vars)
    func = eval(Expr(:->, Expr(:parameters,keys(vars)...),expr))
    return Base.invokelatest(func;vars...)
end

function TTM_simulation(du,u,p,t)
    Tel_idx = findfirst(==("Tel"),p[2]) 
    Tph_idx = findfirst(==("Tph"),p[2]) 
    μ = find_chemicalpotential(p[6],u.x[Tel_idx][1],p[4].DOS,p[5].kB)
    vars = (Tel=u.x[Tel_idx][1],Tph=u.x[Tph_idx][1],t=t,μ=μ,las=p[3],cons=p[5],mp=p[4]) 
    Threads.@threads for i in eachindex(p[2])
        du.x[i][:] .= eval_expr(p[1][p[2][i]],vars)
    end
end

function AthEM_simulation(du,u,p,t)
    println(t)
    Tel_idx = findfirst(==("Tel"),p[2])
    neq_idx = findfirst(==("fneq"),p[2])
    n_idx = findfirst(==("noe"),p[2])

    μ = find_chemicalpotential(u.x[n_idx][1],u.x[Tel_idx][1],p[4].DOS,p[5].kB)

    relax_vars = (Tel=u.x[Tel_idx][1],fneq=u.x[neq_idx][:],n=u.x[n_idx][1],μ=μ,cons=p[5],mp=p[4],n0=p[6],u0=p[7],t=t,las=p[3])
    relax_dis = eval_expr(athem_electronelectronscattering(),relax_vars)

    n_vars = merge(relax_vars,(Symbol("relax_dis")=>relax_dis,))
    du.x[n_idx][:] .= eval_expr(p[1]["noe"],n_vars)

    vars = merge(n_vars,(Symbol("Δn")=>du.x[n_idx][1][1],))
    du.x[Tel_idx][:].=eval_expr(p[1]["Tel"],vars)
    du.x[neq_idx][:].=eval_expr(p[1]["fneq"],vars)
end
   
function run_dynamics(p,u0::ArrayPartition{Float64,Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}})
    prob=ODEProblem(AthEM_simulation,u0,(-75,200.0),p)
    sol = solve(prob,Tsit5(),abstol=1e-4,reltol=1e-4,saveat=1.0,dtmin=0.1)
    return sol
end

function run_dynamics(p,u0::ArrayPartition{Float64, Tuple{Vector{Float64}, Vector{Float64}}})
    prob=ODEProblem(TTM_simulation,u0,(-75,200.0),p)
    sol = solve(prob,Tsit5(),abstol=1e-2,reltol=1e-2,saveat=1.0,dtmin=0.1)
    return sol
end

function main()
    #Build simulation with settings
    sim,mp,las,dim,cons = setup()
    laser=laser_factory(las,dim)
    exprs,key_list = generate_expressions(sim,laser,dim)

    u0 = get_u0(mp.DOS,0.0)
    n0=get_n0(mp.DOS,0.0)

    p=(exprs,key_list,las,mp,cons,n0,u0)
    u = ArrayPartition([get_thermalparticles(mp.μ,1e-16,mp.DOS,cons.kB)],[300.0],zeros(length(mp.egrid)))#should match the order of key_list

    sol=run_dynamics(p,u)
    return sol
end

sol=main()
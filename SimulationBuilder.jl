using ModelingToolkit,DifferentialEquations,Plots,Symbolics,Dierckx,DelimitedFiles,Integrals
using Unitful,BenchmarkTools,ForwardDiff,StaticArrays
using ModelingToolkit: t_nounits as t, D_nounits as D 

include("SymbolicsInterpolation.jl")
include("SimulationVariables.jl")
include("SimulationSetup.jl")
include("Lasers.jl")
include("ElectronTemperature.jl")
include("PhononTemperature.jl")
include("ElectronDistribution.jl")


function setup()
    sim = define_simulation_settings(nlchempot=true,nlelecphon=true,nlelecheat=true,noneqelec=false,elecphonint=true)
    mp = define_material_parameters(extcof=12.7,gamma=4.4315e-22,debye=165,noatoms=59,plasma=2.1357,thermalcond=320.0,
    elecperatom=1,eleceffmass=1.1,dos="DOS/Au_DOS.dat",secmomspecfun=23e-6,elecphon=1.44e-7,ballistic=0.0,cph=0.015)
    las=define_laser_system(lasertype="Gaussian",fwhm=50,fluence=62.42,photon_en=3.1,offset=200)
    laser=laser_factory(las)
    return sim,mp,las,laser
end

#= function setup()
    sim = define_simulation_settings(noneqelec=false,phononheatcapacity=false)
    mp = define_material_parameters(extcof=12.7,gamma=4.4315e-7,debye=165,noatoms=59,plasma=2.1357,thermalcond=320.0,
    elecperatom=1,eleceffmass=1.1,dos="DOS/Au_DOS.dat",secmomspecfun=23e-6,elecphon=1.44e-7,ballistic=0.0,cph=0.015)
    las=define_laser_system(lasertype="Gaussian",fwhm=50,fluence=62.42,photon_en=3.1,offset=200)
    laser=laser_factory(las)
    return sim,mp,las,laser
end  =#

function equation_builder(sim,mp,laser)
    @named Phonon_temp=t_phonon_factory(mp,sim)
    @named Electron_temp = t_electron_factory(mp,sim,laser)
    no_part = get_thermalparticles(0.0,1e-16,mp.DOS,8.617e-5)
    connections=[Electron_temp.Tph ~ Phonon_temp.Tph,
                Electron_temp.Tel ~ Phonon_temp.Tel]
    chempot = (t>=0.0) => (update_chempot!,[Electron_temp.dTel.Tel=>:Tel],
    [Electron_temp.μ=>:μ,Electron_temp.kB=>:kB],[Electron_temp.μ],(mp.DOS,no_part))
    connected = compose(ODESystem(connections,t,name=:connected,defaults=Pair{Num,Any}[Phonon_temp.kB => Electron_temp.kB
    ,Phonon_temp.λ=>Electron_temp.λ,Phonon_temp.hbar=>Electron_temp.hbar,Phonon_temp.μ=>Electron_temp.μ],discrete_events=chempot)
    ,Electron_temp,Phonon_temp)
    connected_simp=structural_simplify(connected)
    return connected_simp,Electron_temp,Phonon_temp
end

#= function equation_builder(sim,mp,laser)
    @named Phonon_temp=t_phonon_factory(mp,sim)
    @named Electron_temp = t_electron_factory(mp,sim,laser)
    connections=[Electron_temp.Tph ~ Phonon_temp.Tph,
                Electron_temp.Tel ~ Phonon_temp.Tel]
    connected = compose(ODESystem(connections,t,name=:connected,defaults=Pair{Num,Any}[Phonon_temp.g => Electron_temp.g]),Electron_temp,Phonon_temp)
    connected_simp=structural_simplify(connected)
    return connected_simp,Electron_temp,Phonon_temp
end =#

function run_dynamics(connected_eq,Tel_eq,Tph_eq,las,mp)
    u0=[Tel_eq.Tel=>300.0,
        Tph_eq.Tph=>300.0]
    p=[Tel_eq.kB => 8.617e-5,
    Tel_eq.μ => 0.0,
    Tel_eq.λ => mp.λ,
    Tel_eq.hbar => 0.6582,
    Tel_eq.FWHM => las.FWHM,
    Tel_eq.Offset => las.Offset,
    Tel_eq.ϵ => mp.ϵ,
    Tel_eq.ϕ => las.Power,
    Tel_eq.R => las.Reflectivity,
    Tph_eq.n => mp.n,
    Tph_eq.θ => mp.θ]
    #= p=[Tel_eq.FWHM => las.FWHM,
    Tel_eq.Offset => las.Offset,
    Tel_eq.ϵ => mp.ϵ,
    Tel_eq.ϕ => las.Power,
    Tel_eq.R => las.Reflectivity,
    Tel_eq.g => mp.g,
    Tel_eq.γ => mp.γ] =#

    tspan=(0.0,500.0)
    prob=ODEProblem(connected_eq,u0,tspan,p)
    sol=solve(prob,Tsit5();abstol=1e-3,reltol=1e-3)
    return sol
end

function main()
    sim,mp,las,laser=setup()
    connected_eq,Tel_eq,Tph_eq=equation_builder(sim,mp,laser)
    sol=run_dynamics(connected_eq,Tel_eq,Tph_eq,las,mp)
    plot(sol)
end

@btime main()
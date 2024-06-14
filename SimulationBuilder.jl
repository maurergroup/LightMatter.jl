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

sim = define_simulation_settings(nlchempot=true,nlelecphon=true,nlelecheat=true,noneqelec=false,elecphonint=true)
mp = define_material_parameters(extcof=12.7,gamma=71,debye=165,noatoms=59,plasma=2.1357,thermalcond=320.0,
elecperatom=1,eleceffmass=1.1,dos="DOS/Au_DOS.dat",secmomspecfun=23e6,elecphon=1.44e-7,ballistic=0.0,cph=0.0)
@named Phonon_temp=t_phonon_factory(mp,sim)
las=define_laser_system(lasertype="Gaussian",fwhm=50,fluence=62.42,photon_en=3.1,offset=200)
laser=laser_factory(las)
@named Electron_temp = t_electron_factory(mp,sim,laser)

connections=[Electron_temp.Tph ~ Phonon_temp.Tph,
            Electron_temp.Tel ~ Phonon_temp.Tel]
connected = compose(ODESystem(connections,t,name=:connected
                ,defaults=Pair{Num,Any}[Electron_temp.kB => Phonon_temp.kB,Phonon_temp.λ=>Electron_temp.λ,Phonon_temp.hbar=>Electron_temp.hbar]),Electron_temp,Phonon_temp)
connected_simp=structural_simplify(connected)

u0=[Electron_temp.Tel=>300.0,
    Phonon_temp.Tph=>300.0]
p=[Phonon_temp.kB => 8.617e-5,
   Electron_temp.μ => 0.0,
   Electron_temp.λ => mp.λ,
   Electron_temp.hbar => 0.6582,
   Electron_temp.FWHM => las.FWHM,
   Electron_temp.Offset => las.Offset,
   Electron_temp.ϵ => mp.ϵ,
   Electron_temp.ϕ => las.Power,
   Electron_temp.R => las.Reflectivity,
   Phonon_temp.n => mp.n,
   Phonon_temp.θ => mp.θ,
   Phonon_temp.g => mp.g]

tspan=(0.0,500.0)
prob=ODEProblem(connected_simp,u0,tspan,p)
sol=solve(prob,Tsit5();abstol=1e-5,reltol=1e-5)
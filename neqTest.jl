using ModelingToolkit,DifferentialEquations,Plots,Symbolics,Dierckx,DelimitedFiles,Integrals
using Unitful,BenchmarkTools,ForwardDiff,StaticArrays,IfElse
using ModelingToolkit: t_nounits as t, D_nounits as D 

include("SymbolicsInterpolation.jl")
include("SimulationVariables.jl")
include("SimulationSetup.jl")
include("Lasers.jl")
include("ElectronTemperature.jl")
include("PhononTemperature.jl")
include("ElectronDistribution.jl")

DOS=generate_DOS("DOS/Au_DOS.dat",59)
ERange=collect(range(-10,10,step=0.1))
las=define_laser_system(:Gaussian,fwhm=50,fluence=62.42,photon_en=3.1)
laser=laser_factory(las)
@named test = athem_factory(DOS,laser,length(ERange))
test_comp=structural_simplify(test)
u0=[neqFD => zeros(201)]
prob=ODEProblem(test_comp,u0,(0,500))
sol=solve(prob,Tsit5();abstol=1e-3,reltol=1e-3)
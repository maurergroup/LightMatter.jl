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

function sim_builder()
    DOS=generate_DOS("DOS/Au_DOS.dat",59)
    ERange=collect(range(-10,-9,step=0.1))
    las=define_laser_system(:Gaussian,fwhm=50,fluence=62.42,photon_en=3.1)
    laser=laser_factory(las)
    @named neq_eq = athem_factory(DOS,laser,length(ERange))
    neq_simp=structural_simplify(neq_eq)
    return neq_simp, neq_eq,ERange
end

function run_dynamics(neq_simp,neq,ERange)
    u0 = [neq_simp.dneqFD.fneq => zeros(length(ERange))]
    p = [neq_simp.μ => 0.0,
        neq_simp.egrid => ERange,
        neq_simp.FWHM => 50.0,
        neq_simp.kB => 8.617e-5,
        neq_simp.hv => 3.1,
        neq_simp.ϵ => 12.7,
        neq_simp.ϕ => 62.41,
        neq_simp.Tel => 300.0,
        neq_simp.R => 0.0]
    
    tspan = (-200.0,-199.0)
    prob=ODEProblem(neq_simp,u0,tspan,p)
    #sol = solve(prob,Tsit5();abstol=1e-3,reltol=1e-3)
    #return prob
end
function main()
   neq_simp,neq,ERange = sim_builder()
   sol = run_dynamics(neq_simp,neq,ERange) 
end

eqs=main()
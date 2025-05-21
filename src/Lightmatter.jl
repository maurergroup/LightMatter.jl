module Lightmatter

using DataInterpolations, DelimitedFiles, Integrals, NonlinearSolve, RecursiveArrayTools, OrdinaryDiffEq, HDF5, LinearAlgebra, GeneralizedGenerated, Unitful, JLD2, ForwardDiff
using SparseConnectivityTracer, ADTypes

export build_Simulation, run_simulation, post_production, Constants, DensityMatrix, build_DensityMatrix, build_Dimension, FE_initialization
export ElectronicTemperature, build_ElectronicTemperature, PhononicTemperature, build_PhononicTemperature, function_builder
export Laser, build_Laser, AthermalElectrons, build_AthermalElectrons, Structure, build_Structure, get_FermiEnergy
export ElectronicDistribution, PhononicDistribution


include("UnitManagement.jl")
using .UnitModule
"""
    __init__()

    Initialises custom unit functions
"""
function __init__()
    Unitful.register(UnitModule)
    Unitful.preferunits(UnitModule.eVm, u"nm", u"fs", u"K")
end

"""
    Lightmatter_units

    A list of the preferred units in the Lightmatter.jl package
"""
global const Lightmatter_units = [u"eV", u"nm", u"fs", u"K",UnitModule.eVm]

include("SimulationTypes.jl")
include("Lasers.jl")
include("ElectronicTemperature.jl")
include("PhononicTemperature.jl")
include("PropertyFunctions.jl")
include("DOS_Geometry.jl")
include("ElectronicDistribution.jl")
include("AthermalElectrons.jl")
include("PhononicDistribution.jl")
include("DensityMatrix.jl")
include("RunDynamics.jl")
include("SystemConstruction.jl")
include("OutputProcessing.jl")
include("AntennaReactor.jl")

end

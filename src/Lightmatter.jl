module Lightmatter

using DataInterpolations, DelimitedFiles, Integrals, NonlinearSolve, RecursiveArrayTools, OrdinaryDiffEq, HDF5, LinearAlgebra, GeneralizedGenerated, Unitful, JLD2

export build_Simulation, run_simulation, post_production, Constants, DensityMatrix, build_DensityMatrix, build_Dimension, FE_initialization
export ElectronicTemperature, build_ElectronicTemperature, PhononicTemperature, build_PhononicTemperature, function_builder
export Laser, build_Laser, AthermalElectrons, build_AthermalElectrons, Structure, build_Structure, get_FermiEnergy
export ElectronicDistribution, PhononicDistribution


Unitful.uconvert(a::Unitful.FreeUnits, b::Union{Real, Array{<:Real}}) = b

"""
    Lightmatter_units

    A list of the units used in Lightmatter.jl: Please convert all units to this 
"""
global const Lightmatter_units = [u"eV", u"nm", u"fs", u"K"] 

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

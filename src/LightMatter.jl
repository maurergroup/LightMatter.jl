module LightMatter

using DelimitedFiles, Integrals, NonlinearSolve, RecursiveArrayTools, OrdinaryDiffEq, HDF5, LinearAlgebra, GeneralizedGenerated, Unitful, JLD2, ForwardDiff
using Bessels, FillArrays, ArrayInterface, Interpolations, DataInterpolations, PreallocationTools, StaticArrays, SimpleNonlinearSolve

export build_Simulation, run_simulation, post_production, Constants, DensityMatrix, build_DensityMatrix, build_Dimension, FE_initialization
export ElectronicTemperature, build_ElectronicTemperature, PhononicTemperature, build_PhononicTemperature, function_builder
export Laser, build_Laser, AthermalElectrons, build_AthermalElectrons, Structure, build_Structure
export ElectronicDistribution, PhononicDistribution, BaseUnits


Unitful.uconvert(a::Unitful.FreeUnits, b::Union{Real, Array{<:Real}}) = b

"""
    LightMatter_units

    A list of the units used in LightMatter.jl: Please convert all units to this 
"""
global const LightMatter_units = [u"eV", u"nm", u"fs", u"K"] 

"""
    spl=DataInterpolations.LinearInterpolation

    A convenience type definition to make type specificity easier throughout the code
"""
global const spl=Interpolations.AbstractExtrapolation

import Base.getindex

include("UnitManagement.jl")

include("SimulationTypes.jl")
getindex(obj::Simulation, x) = obj

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
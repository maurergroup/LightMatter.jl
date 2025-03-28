"""
    Type that all Simulation settings types will be a subtype of
"""
abstract type Simulation end
"""
    The supertype for the different laser types. Each laser type requires a struct which holds the laser parameters.
    Each laser type also has a functor defined for it which returns the equation describing how the laser evolves 
    over time. Lasers are centred on t = 0
"""
abstract type Laser end
"""
    A convenience type definition to make type specificity easier throughout the code
"""
spl=DataInterpolations.LinearInterpolation
"""
    struct Interaction <:Simulation
        ElectronElectron::Bool # Provides whether electron electron interactions are enabled
        ElectronPhonon::Bool   # Provides whether electron phonon interactions are enabled
    end
"""
@kwdef struct Interaction <:Simulation
    ElectronElectron::Bool
    ElectronPhonon::Bool
end
"""
    struct ParameterApproximation <:Simulation
        ElectronPhononCoupling::Bool # Sets the electron-phonon coupling to non-linear from constant
        ElectronHeatCapacity::Bool # Sets the electron heat capacity to non-linear from constant
        PhononHeatCapacity::Bool # Sets the phonon heat capacity to non-linear from constant
        EmbeddingMethod::Bool # Sets the higher order method (e.g. AthEM) to only the top level of the slab
    end

    True represents that the non-linear forms are enabled. 
"""
@kwdef struct ParameterApproximation <:Simulation
    ElectronPhononCoupling::Bool
    ElectronHeatCapacity::Bool
    PhononHeatCapacity::Bool
    EmbeddingMethod::Bool
end
"""
    Booleans for whether a system should be included in the coupled ODE's. There will be an ODE
    for each true statement within this list.
"""
@kwdef struct SystemComponents <:Simulation
    ElectronTemperature::Bool
    PhononTemperature::Bool
    NonEqElectrons::Bool
end
"""
    The overarching struct that holds all Boolean flags for setting up a simulation as well as further
    information such as simulation length, solver settings and output file names etc.
"""
@kwdef struct SimulationSettings <: Simulation
    ParameterApprox::ParameterApproximation
    Interactions::Interaction
    Systems::SystemComponents
    Spatial_DOS::Bool
    DistributionConductivity::Bool
end
"""
    Dimension is the parent type of all dimensionality information
"""
abstract type Dimension <: Simulation end
"""
    A type to define the model is in 0D
"""
@kwdef struct Homogeneous <: Dimension 
    length::Int
    grid::Vector{Float64}
end
"""
    1D struct that holds a grid defined by it's length with pre-determined spacing
"""
@kwdef struct Linear <: Dimension # The 1D model defined by a depth
    grid::Vector{Float64}
    dz::Real
    length::Int
end
"""
    2D struct that holds a grid defined by it's length and radius with pre-determined spacing
"""
@kwdef struct Cylindrical <: Dimension #This assumes circular symmetry around the laser pulse so x & y are defined by a radius
    spatial_grid::Matrix{<:Real}
end
"""
    3D struct that holds a grid defined by the length of each axis with pre-determined spacing
"""
@kwdef struct Cubic{x,y,z} <: Dimension #A fully 3D cuboidal model but current construction has the length of x & y equal
    spatial_grid::Array{<:Real}
end
"""
    This struct holds all parameters and informaion pertaining to the material. Undecided on whether alloys
    should define multiple material parameters or whether it should become an array.
"""
@kwdef struct MaterialParameters <: Simulation #Holds all material parameters
    ϵ::Float64 #Extinction coefficient
    μ::Float64 # centred at 0.0
    FE::Float64 # centred at the difference between top and bottom of valence bands for FLT relaxation
    γ::Float64 #Specific electronic heat capacity
    θ::Float64 #Debye temperature
    n::Float64 #Number of atoms per volume
    κ::Float64 #Room temperature thermal conductivity
    DOS::Vector{spl} #The DOS
    λ::Float64 #Second momentum of spectral function
    g::Float64 # Linear electron-phonon coupling constant
    δb::Float64 # Ballistic length of electrons
    Cph::Float64 #Constant heat capacity for phonons
    egrid::Vector{Float64} # Energy grid to solve neq electrons on
    τ::Float64 #Scalar value for the Fermi Liquid Theory relaxation time
    n0::Vector{Float64}
    τep::Float64
    R::Real # Reflectivity of the sample
    v_g::Vector{<:Real}
end
"""
    Struct that holds constants
"""
struct Constants
    kB::Float64
    hbar::Float64
    me::Float64
end
"""
    Generates the simulation_settings struct with user inputs and defaults or user settings and a dictionary generated from an input
    file built within InputFileControl.jl. Temporary : File Control not fully supported
"""
function define_simulation_settings(;nlelecphon=false,nlelecheat=false,noneqelec=false,elecelecint=false,elecphonint=false,
    nlphonheat=false,electemp=false,phonontemp=false,zDOS=false,embedding=false,distributionconductivity=false)
    
    params=ParameterApproximation(ElectronPhononCoupling=nlelecphon,ElectronHeatCapacity=nlelecheat,
    PhononHeatCapacity=nlphonheat,EmbeddingMethod=embedding)

    interact=Interaction(ElectronElectron=elecelecint,ElectronPhonon=elecphonint)
    
    components=SystemComponents(ElectronTemperature=electemp,PhononTemperature=phonontemp,NonEqElectrons=noneqelec)

    sim_settings=SimulationSettings(ParameterApprox=params,Interactions=interact,Systems=components,Spatial_DOS=zDOS
    ,DistributionConductivity = distributionconductivity)
    
    return sim_settings
end
"""
    Builds the correct slab based on purely user settings or user settings and a dictionary generated from an input
    file built within InputFileControl.jl. Temporary : File Control not fully supported
"""
function define_sim_dimensions(;Dimension::Int64,Lengths=400::Union{Int,Vector{Int}},spacing=1::Union{Real,Vector{Real}})
    
    if Dimension == 0
        return Homogeneous(length=1,grid=[0.0])
    elseif Dimension == 1
        return one_dimension(Lengths,spacing)
    elseif Dimension == 2
        return two_dimension(Lengths,spacing)
    elseif Dimension == 3
        return three_dimension(Lengths,spacing)
    end
end
"""
    Builds 1D slab using the length and spacing
"""
function one_dimension(Length,spacing)
    zs=collect(0:spacing:Length)
    return Linear(grid=zs,dz=spacing,length=length(zs))
end
"""
    Builds 2D slab, the vectors are ordered by [z-axis,radius]
"""
function two_dimension(Length,spacing)
    rs=Vector(0:spacing[2]:Length[2])
    zs=Vector(0:spacing[1]:Length[1])
    slab_grid=Matrix{Tuple{Float64,Float64}}(undef,length(zs),length(rs))
    for (i,r) in enumerate(rs)
        for (j,z) in enumerate(zs)
            slab_grid[j,i]=(r,z)
        end
    end
    return Cylindrical(spatial_grid=slab_grid)
end
"""
    Builds 3D slab, the vectors are ordered by [z-axis,x-axis,y-axis]
"""
function three_dimension(Length,spacing)
    xs=Vector(-Length[2]/2:spacing[2]:Length[2]/2)
    ys=Vector(-Length[3]/2:spacing[3]:Length[3]/2)
    zs=Vector(0:spacing[1]:Length[1])
    slab_grid=Array{Tuple{Float64,Float64,Float64}}(undef,length(ys),length(xs),length(zs))
    for (i,x) in enumerate(xs)
        for (j,y) in enumerate(ys)
            for (k,z) in enumerate(zs)
                slab_grid[j,i,k]=(x,y,z)
            end
        end
    end
    return Cubic(spatial_grid=slab_grid)
end
"""
    Builds the material parameter struct with user parameters or user settings and a dictionary generated from an input
    file built within InputFileControl.jl Temporary : File Control not fully supported
"""
function define_material_parameters(las::Laser,sim::SimulationSettings,dim::Dimension;extcof=0.0,gamma=0.0,debye=0.0
    ,noatoms=0.0,plasma=0.0,thermalcond=0.0,dos="",bulk=true,secmomspecfun=0.0,elecphon=0.0,ballistic=0.0,cph=0.0,τf=18.0,
    folder=Nothing,surfacegeometry="",bulkgeometry="",layer_tolerance=0.1,skip=0,reflectivity=0.0)
    
    fermien=get_FermiEnergy(dos,skip)
    Vbulk = get_unitcellvolume(bulkgeometry)
    if sim.Spatial_DOS == true
        DOS = spatial_DOS(folder,surfacegeometry,dos,Vbulk,noatoms,dim,layer_tolerance,skip)
    elseif sim.Spatial_DOS == false
        if bulk 
            DOS = fill(generate_DOS(dos,1/Vbulk,skip),dim.length)
        else 
            DOS = fill(generate_DOS(dos,noatoms,skip),dim.length)
        end
    end
    tau = 128/(sqrt(3)*pi^2*plasma)
    erange = build_egrid(las.hv)
    n0 = zeros(dim.length)
    for i in eachindex(n0)
        n0[i] = get_thermalparticles(0.0,1e-32,DOS[i],8.617e-5,erange)
    end
    τep = τf*las.hv/8.617e-5/debye

    v_g = get_fermigas_velocity(erange,fermien)

    matpat=MaterialParameters(ϵ=extcof,μ=0.0,γ=gamma,θ=debye,n=noatoms,κ=thermalcond,v_g = v_g,
    DOS=DOS,λ=secmomspecfun,g=elecphon,δb=ballistic,Cph=cph,egrid=erange,τ = tau,FE=fermien,n0=n0,τep=τep,R=reflectivity)
    return matpat
end

function build_egrid(hv)
    egrid = collect(range(-2*hv,2*hv,step=0.01))
    side = 1
    while length(egrid) % 4 != 1
        if side == 1
            push!(egrid,egrid[end]+0.01)
            side = -1
        elseif side == -1
            pushfirst!(egrid,egrid[1]-0.01)
            side = 1
        end
    end
    return egrid
end

function get_fermigas_velocity(egrid,EF)
    eV_to_J(E) = E * 1.602e-19
    ms_to_nmfs(v) = v*1e-6
    v_g = ms_to_nmfs.(sqrt.(2*eV_to_J.(egrid.+EF)./cons.me))
    return v_g
end

const cons=Constants(8.617e-5,0.6582,3.109e-31)
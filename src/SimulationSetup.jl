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
    Booleans for interactions between different ODE systems are held within this struct.
    They are used as flags to determine whether couplings should evaluate to a function
    if true or 0.0 if set to false
"""
@kwdef struct Interaction <:Simulation
    ElectronElectron::Bool
    ElectronPhonon::Bool
end
"""
    Booleans for whether parameters are variable or static. In the case, of the chemical potential
    it is updated via a struct whereas the rest are flags for function calls. If set to true then
    the non-linear/updating form is used.
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
end
"""
    Struct that holds constants
"""
struct Constants
    kB::Float64
    hbar::Float64
end
"""
    Generates the simulation_settings struct with user inputs and defaults or user settings and a dictionary generated from an input
    file built within InputFileControl.jl. Temporary : File Control not fully supported
"""
function define_simulation_settings(;nlelecphon=false,nlelecheat=false,noneqelec=true,elecelecint=true,elecphonint=true,
    phononheatcapacity=true,electemp=true,phonontemp=true,zDOS=false,embedding=false)
    
    if noneqelec==false
        elecelecint=false
    end

    params=ParameterApproximation(ElectronPhononCoupling=nlelecphon,ElectronHeatCapacity=nlelecheat,
    PhononHeatCapacity=phononheatcapacity,EmbeddingMethod=embedding)

    interact=Interaction(ElectronElectron=elecelecint,ElectronPhonon=elecphonint)
    
    components=SystemComponents(ElectronTemperature=electemp,PhononTemperature=phonontemp,NonEqElectrons=noneqelec)

    sim_settings=SimulationSettings(ParameterApprox=params,Interactions=interact,Systems=components,Spatial_DOS=zDOS)
    
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
    ,noatoms=0.0,plasma=0.0,thermalcond=0.0,dos="",secmomspecfun=0.0,elecphon=0.0,ballistic=0.0,cph=0.0,τf=18.0,
    folder=Nothing,surfacegeometry="",bulkgeometry="",layer_tolerance=0.1,skip=0,reflectivity=0.0)
    
    fermien=get_FermiEnergy(dos,skip)
    Vbulk = get_unitcellvolume(bulkgeometry)
    if sim.Spatial_DOS == true
        DOS = spatial_DOS(folder,surfacegeometry,dos,Vbulk,noatoms,dim,layer_tolerance,skip)
    elseif sim.Spatial_DOS == false
        if typeof(dim) == Homogeneous
            DOS = [generate_DOS(dos,Vbulk,skip)]
        else
            DOS = fill(generate_DOS(dos,Vbulk,skip),dim.length)
        end
    end
    tau = 128/(sqrt(3)*pi^2*plasma)
    erange = build_egrid(las.hv)#
    n0 = zeros(dim.length)
    for i in eachindex(n0)
        n0[i] = get_thermalparticles(0.0,1e-32,DOS[i],8.617e-5,erange)
    end
    τep = τf*las.hv/8.617e-5/debye

    matpat=MaterialParameters(ϵ=extcof,μ=0.0,γ=gamma,θ=debye,n=noatoms,κ=thermalcond,
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

function grid_builder(l,Espan)
    gh,weights = gausshermite(l)
    return ((gh .- minimum(gh)) ./ (maximum(gh)/(Espan/2))) .- (Espan/2) 
end

const cons=Constants(8.617e-5,0.6582)
"""
    Type that all abstract types will derive from
"""
abstract type Simulation end #Overarching type that holds all simulation settings
abstract type Laser <: Simulation end
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
    dz::Int
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
    ne::Float64 #Number of electrons per atom
    effmass::Float64 #Effective mass of conduction electrons
    DOS::spl #The DOS
    λ::Float64 #Second momentum of spectral function
    g::Float64 # Linear electron-phonon coupling constant
    Ballistic::Float64 # Ballistic length of electrons
    Cph::Float64 #Constant heat capacity for phonons
    egrid::Vector{Float64} # Energy grid to solve neq electrons on
    τ::Float64 #Scalar value for the Fermi Liquid Theory relaxation time
    u0::Float64
    n0::Float64
    τep::Float64
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
    phononheatcapacity=true,electemp=true,phonontemp=true)
    
    if noneqelec==false
        elecelecint=false
    end

    params=ParameterApproximation(ElectronPhononCoupling=nlelecphon,ElectronHeatCapacity=nlelecheat,
    PhononHeatCapacity=phononheatcapacity)

    interact=Interaction(ElectronElectron=elecelecint,ElectronPhonon=elecphonint)
    
    components=SystemComponents(ElectronTemperature=electemp,PhononTemperature=phonontemp,NonEqElectrons=noneqelec)

    sim_settings=SimulationSettings(ParameterApprox=params,Interactions=interact,Systems=components)
    
    return sim_settings
end
"""
    Builds the correct slab based on purely user settings or user settings and a dictionary generated from an input
    file built within InputFileControl.jl. Temporary : File Control not fully supported
"""
function define_sim_dimensions(;Dimension::Int64,Lengths=400::Union{Int,Vector{Int}},spacing=1::Union{Int,Vector{Int}})
    
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
function define_material_parameters(las::Laser;extcof=0.0,gamma=0.0,debye=0.0,noatoms=0.0,plasma=0.0,thermalcond=0.0,elecperatom=0.0,eleceffmass=0.0
    ,dos="DOS/Au_DOS.dat",secmomspecfun=0.0,elecphon=0.0,ballistic=0.0,cph=0.0,τf=18.0)
    
    fermien=get_FermiEnergy(dos)
    DOS = generate_DOS(dos,noatoms)
    tau = 128/(sqrt(3)*pi^2*plasma)
    erange = grid_builder(0.0,-3*las.hv,3*las.hv,0.0005, 0.0002) 
    u0 = get_u0(DOS,0.0,fermien)
    n0 = get_n0(DOS,0.0,fermien)
    τep = τf*las.hv/8.617e-5/debye

    matpat=MaterialParameters(ϵ=extcof,μ=0.0,γ=gamma,θ=debye,n=noatoms,κ=thermalcond,ne=elecperatom,effmass=eleceffmass,
    DOS=DOS,λ=secmomspecfun,g=elecphon,Ballistic=ballistic,Cph=cph,egrid=erange,τ = tau,FE=fermien,u0=u0,n0=n0,τep=τep)

    return matpat
end

function grid_builder(mid, left, right, step, α)
    v = Vector{typeof(mid + α*step)}()
    let s = step, a = mid - s
        while a ≥ left
            push!(v, a)
            s += α
            a -= s
        end
    end
    reverse!(v)
    push!(v, mid)
    let s = step, a = mid + s
        while a ≤ right
            push!(v, a)
            s += α
            a += s
        end
    end
    v
end
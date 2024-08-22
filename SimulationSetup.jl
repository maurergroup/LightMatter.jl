"""
    Type that all abstract types will derive from
"""
abstract type Simulation end #Overarching type that holds all simulation settings
abstract type Laser <: Simulation end # Holds all Lasers
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
    ChemicalPotential::Bool
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
    Output_File_Name::String
    Simulation_End_Time::Real
end
"""
    Dimension is the parent type of all dimensionality information
"""
abstract type Dimension <: Simulation end
"""
    A type to define the model is in 0D
"""
struct Homogenous <: Dimension end
"""
    1D struct that holds a grid defined by it's length with pre-determined spacing
"""
@kwdef struct Linear{L} <: Dimension # The 1D model defined by a depth
    spatial_grid::SVector{L,<:Real}
end
"""
    2D struct that holds a grid defined by it's length and radius with pre-determined spacing
"""
@kwdef struct Cylindrical{r,L} <: Dimension #This assumes circular symmetry around the laser pulse so x & y are defined by a radius
    spatial_grid::SMatrix{r,L,Real}
end
"""
    3D struct that holds a grid defined by the length of each axis with pre-determined spacing
"""
@kwdef struct Cubic{x,y,z} <: Dimension #A fully 3D cuboidal model but current construction has the length of x & y equal
    spatial_grid::SArray{Tuple{y,x,z},Real}
end
"""
    This struct holds all parameters and informaion pertaining to the material. Undecided on whether alloys
    should define multiple material parameters or whether it should become an array.
"""
@kwdef struct MaterialParameters <: Simulation #Holds all material parameters
    ϵ::Real #Extinction coefficient
    μ::Real # centred at 0.0
    FE::Real # centred at the difference between top and bottom of valence bands for FLT relaxation
    γ::Real #Specific electronic heat capacity
    θ::Real #Debye temperature
    n::Real #Number of atoms per volume
    κ::Real #Room temperature thermal conductivity
    ne::Real #Number of electrons per atom
    effmass::Real #Effective mass of conduction electrons
    DOS::Spline1D #The DOS
    λ::Real #Second momentum of spectral function
    g::Real # Linear electron-phonon coupling constant
    Ballistic::Real # Ballistic length of electrons
    Cph::Real #Constant heat capacity for phonons
    egrid::Vector{Float64} # Energy grid to solve neq electrons on
    τ::Real #Scalar value for the Fermi Liquid Theory relaxation time
end
"""
    Struct that holds constants
"""
struct Constants
    kB::Real
    hbar::Real
end
"""
    Generates the simulation_settings struct with user inputs and defaults or user settings and a dictionary generated from an input
    file built within InputFileControl.jl. Temporary : File Control not fully supported
"""
function define_simulation_settings(;nlchempot=false,nlelecphon=false,nlelecheat=false,noneqelec=true,elecelecint=true,elecphonint=true,
    phononheatcapacity=true,electemp=true,phonontemp=true,output="./Default_file.jld2",simendtime=1000)
    
    if noneqelec==false
        elecelecint=false
    end

    params=ParameterApproximation(ChemicalPotential=nlchempot,ElectronPhononCoupling=nlelecphon,ElectronHeatCapacity=nlelecheat,
    PhononHeatCapacity=phononheatcapacity)

    interact=Interaction(ElectronElectron=elecelecint,ElectronPhonon=elecphonint)
    
    components=SystemComponents(ElectronTemperature=electemp,PhononTemperature=phonontemp,NonEqElectrons=noneqelec)

    sim_settings=SimulationSettings(ParameterApprox=params,Interactions=interact,Systems=components,Output_File_Name=output,
    Simulation_End_Time=simendtime)
    
    return sim_settings
end

function define_simulation_settings(dict::Dict;nlchempot=dict["ChemPot"],nlelecphon=dict["ElecPhonCoup"],nlelecheat=dict["ElecHeatCapac"],
    noneqelec=dict["NoneqElec"],elecelecint=dict["Elec-Elec"],elecphonint=dict["Elec-Phon"],output=dict["Output"],electemp=dict["Tel"],
    phonontemp=dict["Tph"],simendtime=dict["SimTime"],phononheatcapacity=dict["PhononHeatCapac"])
    
    if noneqelec==false
        if elecelecint==true
            throw(ErrorException("Electron-Electron Interactions should only 
            be true if there are non-equilibrium electrons to interact with. 
            Set elecelecint=false"))
        end
    end

    nl=ParameterApproximation(ChemicalPotential=nlchempot,ElectronPhononCoupling=nlelecphon,ElectronHeatCapacity=nlelecheat,
    PhononHeatCapacity=phononheatcapacity)

    interact=Interaction(ElectronElectron=elecelecint,ElectronPhonon=elecphonint)
    
    components=SystemComponents(ElectronTemperature=electemp,PhononTemperature=phonontemp,NonEqElectrons=noneqelec)

    sim_settings=SimulationSettings(ParameterApprox=nl,Interactions=interact,Systems=components,Output_File_Name=output,
    Simulation_End_Time=simendtime)
    
    return sim_settings
end
"""
    Builds the correct slab based on purely user settings or user settings and a dictionary generated from an input
    file built within InputFileControl.jl. Temporary : File Control not fully supported
"""
function define_sim_dimensions(;Dimension::Int64,Lengths=400::Union{Real,Vector{<:Real}},spacing=1::Union{Real,Vector{<:Real}})
    
    if Dimension == 0
        return Homogenous()
    elseif Dimension == 1
        return one_dimension(Lengths,spacing)
    elseif Dimension == 2
        return two_dimension(Lengths,spacing)
    elseif Dimension == 3
        return three_dimension(Lengths,spacing)
    end
end

function define_sim_dimensions(dict::Dict;Dimension=dict.Dimension::Int64,Lengths=dict.Lengths::Union{Real,Vector{<:Real}}
    ,spacing=dict.Spacing::Union{Real,Vector{<:Real}})

    Dims == length(Lengths) == length(Spacing) || error("The number of dimensions specified doesn't equal the lengths of
    one of the input settongs.")

    if dict.Dimension == 0
        return Homogenous
    elseif dict.Dimension == 1
        return one_dimension(Lengths,spacing)
    elseif dict.Dimension == 2
        return two_dimension(Lengths,spacing)
    elseif dict.Dimension == 3
        return three_dimension(Lengths,spacing)
    end
end
"""
    Builds 1D slab using the length and spacing
"""
function one_dimension(Length,spacing)
    zs=0:spacing:Length
    zs=collect(zs)
    slab = SVector{length(zs),typeof(zs[1])}(zs)
    return Linear{length(slab)}(spatial_grid=slab)
end
"""
    Builds 2D slab, the vectors are ordered by [z-axis,radius]
"""
function two_dimension(Length,spacing)
    rs=SVector{length(0:spacing[2]:Length[2])}(0:spacing[2]:Length[2])
    zs=SVector{length(0:spacing[1]:Length[1])}(0:spacing[1]:Length[1])
    slab_grid=Matrix{Tuple{Real,Real}}(undef,length(zs),length(rs))
    for (i,r) in enumerate(rs)
        for (j,z) in enumerate(zs)
            slab_grid[j,i]=(r,z)
        end
    end
    return Cylindrical{length(rs),length(zs)}(spatial_grid=slab_grid)
end
"""
    Builds 3D slab, the vectors are ordered by [z-axis,x-axis,y-axis]
"""
function three_dimension(Length,spacing)
    xs=SVector{length(-Length[2]/2:spacing[2]:Length[2]/2)}(-Length[2]/2:spacing[2]:Length[2]/2)
    ys=SVector{length(-Length[3]/2:spacing[3]:Length[3]/2)}(-Length[3]/2:spacing[3]:Length[3]/2)
    zs=SVector{length(0:spacing[1]:Length[1])}(0:spacing[1]:Length[1])
    slab_grid=Array{Tuple{Real,Real,Real}}(undef,length(ys),length(xs),length(zs))
    for (i,x) in enumerate(xs)
        for (j,y) in enumerate(ys)
            for (k,z) in enumerate(zs)
                slab_grid[j,i,k]=(x,y,z)
            end
        end
    end
    return Cubic{length(ys),length(xs),length(zs)}(spatial_grid=slab_grid)
end
"""
    Builds the material parameter struct with user parameters or user settings and a dictionary generated from an input
    file built within InputFileControl.jl Temporary : File Control not fully supported
"""
function define_material_parameters(las::Laser;extcof=0.0,gamma=0.0,debye=0.0,noatoms=0.0,plasma=0.0,thermalcond=0.0,elecperatom=0.0,eleceffmass=0.0
    ,dos="DOS/Au_DOS.dat",secmomspecfun=0.0,elecphon=0.0,ballistic=0.0,cph=0.0)
    
    fermien=get_FermiEnergy(dos)
    DOS = generate_DOS(dos,noatoms)
    tau = 0.546#128/(sqrt(3)*pi^2*plasma)
    erange = grid_builder(0.0,-3*las.hv,3*las.hv,0.0005, 0.0002) 
    matpat=MaterialParameters(ϵ=extcof,μ=0.0,γ=gamma,θ=debye,n=noatoms,κ=thermalcond,ne=elecperatom,effmass=eleceffmass,
    DOS=DOS,λ=secmomspecfun,g=elecphon,Ballistic=ballistic,Cph=cph,egrid=erange,τ = tau,FE=fermien)

    return matpat
end

function define_material_parameters(dict::Dict;extcof=dict["ExtCof"],gamma=dict["Gamma"],debye=dict["Debye"],noatoms=dict["AtomDens"],
    plasma=dict["Plasma"],thermalcond=dict["RTKappa"],elecperatom=dict["ne"],eleceffmass=dict["EffMass"],dos=dict["DOS"],
    secmomspecfun=Dict["SpectralFunc"],elecphon=Dict["g"],ballistic=dict["BallisticLength"],cph=dict["Cph"])

    fermien=0.0
    DOS = generate_DOS(dos,fermien)

    matpat=MaterialParameters(ϵ=extcof,FE=fermien,γ=gamma,θ=debye,n=noatoms,ω=plasma,κ=thermalcond,ne=elecperatom,effmass=eleceffmass,
    DOS=DOS,λ=secmomspecfun,g=elecphon,Ballistic=ballistic,Cph,cph)

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
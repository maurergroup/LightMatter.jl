module SimulationSetup
using StaticArrays #For the Dimensionality
abstract type Simulation end #Overarching type that holds all simulation settings

@kwdef struct Interaction <:Simulation #Holds Booleans for how different systems interact
    Nonequilibrium_electrons::Bool
    ElectronElectron::Bool
    ElectronPhonon::Bool
end

@kwdef struct ParameterApproximation <:Simulation #Holds Booleans for the parameter approximations
    ChemicalPotential::Bool
    ElectronPhonon_Coupling::Bool
    ElectronHeatCapacity::Bool
end

@kwdef struct SystemComponents <:Simulation #Which ODE's you want to solve
    ElectronTemperature::Bool
    PhononTemperature::Bool
    NonEqElectrons::Bool
end

@kwdef struct Simulation_Settings <: Simulation #Struct containing all the settings for the simulation
    ParameterApprox::ParameterApproximation
    Interactions::Interaction
    Systems::SystemComponents
    Output_File_Name::String
    Simulation_End_Time::Float64
end

abstract type Dimension end #Overarching type for the dimensionality of the system

struct Homogenous <: Dimension end #A type to define 0D

@kwdef struct Linear{L} <: Dimension # The 1D model defined by a depth
    spatial_grid::SVector{L,Real}
end

@kwdef struct Cylindrical{r,L} <: Dimension #This assumes circular symmetry around the laser pulse so x & y are defined by a radius
    spatial_grid::SMatrix{r,L,Real}
end

@kwdef struct Cubic{x,y,z} <: Dimension #A fully 3D cuboidal model but current construction has the length of x & y equal
    spatial_grid::SArray{Tuple{y,x,z},Real}
end

@kwdef struct MaterialParameters <: Simulation #Holds all material parameters
    ϵ::Float64 #Extinction Coefficient
    FE::Float64 #Fermi Energy
    γ::Float64 #Specific electronic heat capacity
    θ::Float64 #Debye Temperature
    n::Float64 #Number of atoms per volume
    ω::Float64 #Plasma frequency
    κ::Float64 #Room temperature thermal conductivity
    L::Float64 #Length of slab
    dz::Float64 #Distance between points in slab
    ne::Float64 #Number of electrons per atom
    effmass::Float64 #Effective mass of conduction electrons
    DOS::Spline1D #File location of the DOS data
    λ::Float64 #Second momentum of spectral function
    g::Float64 # Linear Electron-phonon coupling constant
    ballistic::Real # Ballistic length of electrons
end

"""
    Generates the simulation_settings struct with user inputs and defaults or user settings and a dictionary generated from an input
    file built within InputFileControl.jl
"""
function define_simulation_settings(;nlchempot=false,nlelecphon=false,
    nlelecheat=false,noneqelec=true,elecelecint=true,elecphonint=true,
    output="./Default_file.jld2",dim=0,simendtime=1000)
    
    if noneqelec==false
        if elecelecint==true
            throw(ErrorException("Electron-Electron Interactions should only 
            be true if there are non-equilibrium electrons to interact with. 
            Set elecelecint=false"))
        end
    end

    nl=Nonlinear(ChemicalPotential=nlchempot,ElectronPhononCoupling=nlelecphon,
    ElectronHeatCapacity=nlelecheat)

    interact=Interactions(NonequilibriumElectrons=noneqelec,
    ElectronElectron=elecelecint,ElectronPhonon=elecphonint)

    sim_settings=simulation_settings(nonlinear=nl,interactions=interact,
    Dimension=dim,Output_file_name=output,Simulation_End_Time=simendtime)
    return sim_settings
end

function define_simulation_settings(dict::Dict;nlchempot=dict["ChemPot"],nlelecphon=dict["ElecPhonCoup"],
    nlelecheat=dict["ElecHeatCapac"],noneqelec=dict["NoneqElec"],elecelecint=dict["Elec-Elec"],
    elecphonint=dict["Elec-Phon"],output=dict["Output"],dim=dict["Dimension"],simendtime=dict["SimTime"])

    if noneqelec==false
        if elecelecint==true
            throw(ErrorException("Electron-Electron Interactions should only 
            be true if there are non-equilibrium electrons to interact with. 
            Set elecelecint=false"))
        end
    end

    nl=Nonlinear(ChemicalPotential=nlchempot,ElectronPhononCoupling=nlelecphon,
    ElectronHeatCapacity=nlelecheat)

    interact=Interactions(NonequilibriumElectrons=noneqelec,
    ElectronElectron=elecelecint,ElectronPhonon=elecphonint)

    sim_settings=simulation_settings(nonlinear=nl,interactions=interact,
    Dimension=dim,Output_file_name=output,Simulation_End_Time=simendtime)
    return sim_settings
end
"""
    Builds the correct slab based on purely user settings or user settings and a dictionary generated from an input
    file built within InputFileControl.jl - all spacings are equal currently
"""
function define_sim_dimensions(;Dimension::Int64,Length=400::Integer,radius=200::Integer,
    xwidth=400::Integer,ywidth=400::Integer,spacing=1::Real)
    if Dimension == 0
        return Homogenous()
    elseif Dimension == 1
        return one_dimension(Length,spacing)
    elseif Dimension == 2
        return two_dimension(Length,spacing,radius)
    elseif Dimension == 3
        return three_dimension(length,xwidth,ywidth,spacing)
    end
end

function define_sim_dimensions(dict::Dict;Dimension=dict.Dimension::Int64,Length=dict.Length::Integer,radius=dict.Radius::Integer,
    xwidth=dict.XWidth::Integer,ywidth=dict.YWidth::Integer,spacing=dict.Spacing::Real)
    dict_overwrite(dict,Dimension=Dimension,Length=Length,radius=radius,xwidth=xwidth,ywidth=ywidth,spacing=spacing)
    if dict.Dimension == 0
        return Homogenous
    elseif dict.Dimension == 1
        return one_dimension(Length,spacing)
    elseif dict.Dimension == 2
        return two_dimension(Length,spacing,radius)
    elseif dict.Dimension == 3
        return three_dimension(length,xwidth,ywidth,spacing)
    end
end
"""
    Builds 1D slab
"""
function one_dimension(Length,spacing)
    zs=0:spacing:Length
    zs=collect(zs)
    slab = SVector{length(zs),typeof(zs[1])}(zs)
    return Linear{length(slab)}(spatial_grid=slab)
end
"""
    Builds 2D slab
"""
function two_dimension(Length,spacing,radius)
    rs=SVector{length(0:spacing:radius)}(0:spacing:radius)
    zs=SVector{length(0:spacing:Length)}(0:spacing:Length)
    slab_grid=Matrix{Tuple{Real,Real}}(undef,length(zs),length(rs))
    for (i,r) in enumerate(rs)
        for (j,z) in enumerate(zs)
            slab_grid[j,i]=(r,z)
        end
    end
    return Cylindrical{length(rs),length(zs)}(spatial_grid=slab_grid)
end
"""
    Builds 3D slab
"""
function three_dimension(Length,xwidth,ywidth,spacing)
    xs=SVector{length(-xwidth/2:spacing:xwidth/2)}(-xwidth/2:spacing:xwidth/2)
    ys=SVector{length(-ywidth/2:spacing:ywidth/2)}(-ywidth/2:spacing:ywidth/2)
    zs=SVector{length(0:spacing:Length)}(0:spacing:Length)
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
    file built within InputFileControl.jl
"""
function define_material_parameters(;extcof,gamma,debye,noatoms,
    plasma,thermalcond,elecperatom,eleceffmass,dos,secmomspecfun,elecphon)

    fermien = get_FermiEnergy(dos)
    DOS = generate_DOS(dos,fermien)

    matpat=MaterialParameters(ϵ=extcof,FE=fermien,γ=gamma
    ,θ=debye,n=noatoms,ω=plasma,κ=thermalcond,ne=elecperatom,
    effmass=eleceffmass,DOS=DOS,λ=secmomspecfun,g=elecphon)

    return matpat
end

function define_material_parameters(dict::Dict;extcof=dict["ExtCof"],gamma=dict["Gamma"],
    debye=dict["Debye"],noatoms=dict["AtomDens"],plasma=dict["Plasma"],
    thermalcond=dict["RTKappa"],elecperatom=dict["ne"],eleceffmass=dict["EffMass"],
    dos=dict["DOS"],secmomspecfun=Dict["SpectralFunc"],elecphon=Dict["g"])

    fermien = get_FermiEnergy(dos)
    DOS = generate_DOS(dos,fermien)

    matpat=MaterialParameters(ϵ=extcof,FE=fermien,γ=gamma
    ,θ=debye,n=noatoms,ω=plasma,κ=thermalcond,ne=elecperatom,effmass=eleceffmass,
    DOS=DOS,λ=secmomspecfun,g=elecphon)

    return matpat
end

export Simulation,MaterialParameters,Dimension,Simulation_Settings
export define_simulation_settings,define_simulation_dimensions,define_material_parameters

end
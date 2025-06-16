"""
    SimulationTypes
    
    Parent type of all subtypes in Lightmatter.jl
"""
abstract type SimulationTypes end

"""
    Laser <: SimulationTypes
        envelope::Symbol = :Gaussian # Currently implemented are :Gaussian, :HyperbolicSecant, :Lorentzian and :Rectangular
        Transport::Symbol # :optical, :ballistic and :combined are the options for how the laser decays into a slab

        FWHM::Number # The Full-Width Half-Maximum of the laser, for rectnagular half the length
        ϕ::Number # The unabsorbed fluence of the laser
        hv::Number # The photon frequency of the laser
        ϵ::Union{Number,Vector{<:Number}} # The inverse of the absorption coefficient
        R::Number # The reflectivity of the sample
        δb::Union{Number,Vector{<:Number}} # The ballistic length of electrons
    end
    
    Struct that contains all laser parameters and any material parameters that affect laser absorption
"""
@kwdef struct Laser <: SimulationTypes
    envelope::Symbol = :Gaussian # Currently implemented are :Gaussian, :HyperbolicSecant, :Lorentzian and :Rectangular
    Transport::Symbol # :optical, :ballistic and :combined are the options for how the laser decays into a slab

    FWHM::Number # The Full-Width Half-Maximum of the laser, for rectnagular half the length
    ϕ::Number # The unabsorbed fluence of the laser
    n::Union{Number, Vector{<:Number}, Matrix{<:Number}, Vector{<:Matrix{Number}}}# The real part of the refractive index of the material
    hv::Union{Number, Matrix{<:Number}} # The photon frequency of the laser
    ϵ::Union{Number, Vector{<:Number}, Vector{<:Vector{Number}}} # The inverse of the absorption coefficient
    R::Number # The reflectivity of the sample
    δb::Union{Number, Vector{<:Number}, Vector{<:Vector{Number}}} # The ballistic length of electrons
end
"""
    build_laser(;envelope=:Gaussian, FWHM=10.0, ϕ=10.0, hv=5.0, Transport=:optical, ϵ=1.0, R=0.0, δb=1.0)

    Outer constructor function to assemble the Laser struct. Can handle unit conversions if the user provides a 
    unitful quantity to the laser. 
    Defaults allow any unneccessary parameters for users simulation to be ignored.

    # Arguments
    - 'envelope': Symbol representing the shape of the envelope, :Gaussian, :HyperbolicSecant, :Lorentzian, :Rectangular
    - 'FWHM': unit = fs:  Full-Width Half-Maximum of the laser pulse or half the duration of the Rectangular laser ± 0.0
    - 'ϕ': unit = eV/nm²: The fluence of the laser 
    - 'hv': unit = eV: The photon energy of the laser
    - 'Transport': The method of spatial transport of the laser, :optical, :ballistic, :combined
    - 'ϵ': unit = nm: The penetration depth of the material (1/α), for the :Optical & :Combined transport
    - 'R': unit = unitless:The reflectivity of the sample surface, leave at 0.0 if your provided fluence is the absorbed fluence
    - 'δb': unit = nm:The ballistic length of electrons, for :ballistic & :combined transport

    # Returns
    - The Laser struct with the user settings and neccessary values converted to the correct units
"""
function build_Laser(;envelope::Symbol = :Gaussian, FWHM::Number = NaN, ϕ::Number = NaN, hv::Union{Number, Matrix{<:Number}} = NaN, Transport::Symbol = :optical,
                      ϵ::Union{Number, Vector{<:Number}, Vector{<:Vector{Number}}} = NaN, R::Number = 0.0, δb::Union{Number, Vector{<:Number}, Vector{<:Vector{Number}}} = NaN,
                      n::Union{Number, Vector{<:Number}, Vector{<:Vector{Number}}} = NaN)
    
    FWHM = convert_units(u"fs", FWHM)
    Power = convert_units(u"eV/nm^2", ϕ)
    if hv isa AbstractMatrix
        hv[:,1] = convert_units(u"eV", hv[:,1])
    else
        hv = convert_units(u"eV", hv)
    end
    ϵ = convert_units(u"nm", ϵ)
    δb = convert_units(u"nm", δb)
    return Laser(envelope=envelope, FWHM=FWHM, ϕ=Power, hv=hv, Transport=Transport, ϵ=ϵ, R=R, δb=δb, n=n)
end
"""
    Dimension <: SimulationTypes
        length::Int # The length of the grid, not the depth of the slab
        grid::AbstractArray{<:Number} # The grid the simulation is solved over
        spacing::Union{Number, Vector{<:Number}} #The spacing between grid points
        InterfaceHeight::Union{Number, Vector{<:Number}} # Height sorted list of the interfaces between materials
    end

    Struct that contains all information regarding the spatial grid that the simulation is performed on.
"""
@kwdef struct Dimension <: SimulationTypes
    length::Int # The length of the grid, not the depth of the slab
    grid::AbstractArray{<:Number} # The grid the simulation is solved over
    spacing::Union{Number, Vector{<:Number}} #The spacing between grid points
    InterfaceHeight::Union{Number, Vector{<:Number}} # Height sorted list of the interfaces between materials
end
"""
    build_Dimension(grid=[0.0]::AbstractArray{<:Number}, cutoff=0.0::Union{Vector{<:Number},Number})

    Outer constructor function to assemble the Dimension struct. The user provides an evenly spaced grid 
    and sorted list of interface heights for antenna-reactor complexes. The user must ensure the length of 
    cutoff = Elemental_System - 1. No unit conversion is performed when assembling this struct.
    Defaults allow any unneccessary parameters for users simulation to be ignored.

    # Arguments
    - 'grid': unit = nm: Vector representing spatial grid. If [0.0] then homogenous (0D) calculation is performed
    - 'cutoff': unit = nm:  Sorted list of all interface heights. Only used when Elemental_System > 1. 

    # Returns
    - The Dimension struct with the users grid and interface heights
"""
function build_Dimension(grid::AbstractArray{<:Number} = [NaN], cutoff::Union{Vector{<:Number},Number} = NaN)
    L = length(grid)
    grid = convert_units(u"nm", grid)
    if L > 1
        spacing = grid[2]-grid[1]
    else
        spacing = 1.0
    end
    return Dimension(length=L, grid=grid, spacing=spacing, InterfaceHeight=cutoff)
end
"""
    Fields <: SimulationTypes
        electric::Expr # Expression for the magnitude of the electric field in the simulation
        magnetic::Expr # Expression for the magnitude of the magnetic field in the simulation
    end

    Struct that contains all information regarding electromagnetic fields in the simulation.
"""
struct Fields <: SimulationTypes
    electric::Union{Expr,Number}
    magnetic::Union{Expr,Number}
end
"""
    TotalFields <: SimulationTypes
        laser::Fields # The laser fields in the simulation
        external::Fields # The external fields in the simulation, e.g. a magnetic field
    end

    Struct that contains both the field generated by the laser and any external fields.
"""
struct TotalFields <: SimulationTypes
    laser::Fields # The laser fields in the simulation
    external::Fields # The external fields in the simulation, e.g. a magnetic field
end
"""
    Structure <: SimulationTypes
        Spatial_DOS::Bool # Whether to vary the DOS with height - if so the DOS becomes a vector

        Elemental_System::Int # The number of elemental systems, if > 1 then each constant and vector
                            # of material parameters needs to become a vector of length=Elemental_System

        DOS::Union{spl,Vector{spl},Vector{Vector{spl}}} # The density of states of the simulation
        egrid::Vector{<:Number} # An energy grid for electronic or phononic distributions to be solved on

        dimension::Union{Dimension} # A struct holding all spatial grid structure (0D or 1D)
         bandstructure::Union{Vector{<:DataInterpolations.AkimaInterpolation}, Vector{<:Vector{DataInterpolations.AkimaInterpolation}},Nothing} 
                        # The band structure of the simulation both in terms of k->E and E-> k
    end

    Struct that contains any spatial information including the DOS, the spatial grid to solve the simulation on and
    the elemental composition of the simulation (e.g. an antenna-reactor system would contain two elemental systems)
"""
@kwdef struct Structure <: SimulationTypes
    Spatial_DOS::Bool # Whether to vary the DOS with height - if so the DOS becomes a vector

    Elemental_System::Int # The number of elemental systems, if > 1 then each constant and vector
                          # of material parameters needs to become a vector of length=Elemental_System

    DOS::Union{spl, Vector{spl}} # The density of states of the simulation
    bandstructure::Union{Vector{<:DataInterpolations.AkimaInterpolation}, Vector{<:Vector{DataInterpolations.AkimaInterpolation}},Nothing} # The band structure of the simulation
    egrid::Vector{<:Number} # An energy grid for electronic or phononic distributions to be solved on

    dimension::Union{Dimension} # A struct holding all spatial grid structure (0D or 1D)
    fields::TotalFields # Any laser and external fields in the simulation
end
"""
    build_Structure(; las::Laser=build_Laser(), Spatial_DOS::Bool = false, Elemental_System::Int = 1, dimension::Dimension = build_Dimension(),
                    bulk_DOS::Union{String,Vector{String},Nothing} = nothing, DOS_folder::Union{String,Vector{String},Nothing} = nothing, 
                    bulk_geometry::Union{String,Vector{String},Nothing} = nothing, slab_geometry::Union{String,Vector{String},Nothing} = nothing, 
                    atomic_layer_tolerance::Union{Number,Vector{Number}} = 0.1, DOS::Union{spl,Vector{spl},Nothing} = nothing, 
                    egrid::Union{Vector{<:Number},Nothing} = nothing)

    Outer constructor function to assemble the Structure struct. No unit conversion is performed.
    All DOS files must be in the format |energy (eV), states (eV⁻¹Vᵤ⁻¹)|. Comment lines (#) are ignored
    Defaults allow any unneccessary parameters for users simulation to be ignored.

    # Arguments
    - 'las': Laser struct, provide if not providing a pre-made energy grid
    - 'Spatial_DOS': Bool for determening whether the DOS is spatially resolved or bulk
    - 'Elemental_System': Number of different crystal systems in the structure
    - 'dimension': Dimension struct, provide if not wanting a 0D calculation
    - 'bulk_DOS': File location of the bulk DOS file
    - 'DOS_folder': Location of a folder containing atom projected DOS. These must be in units of (eV⁻¹atom⁻¹) and be .dat files
    - 'bulk_geometry': File location of the bulk DOS' geometry.in file
    - 'slab_geometry': File location of the geometry.in to create the atom projected DOS' found in DOS_folder
    - 'atomic_layer_tolerance': Minimum height spacing in Å for two atoms to be apart to be considered different layers
    - 'DOS': A user made DOS or vector of DOS' for the simulation if not wanting to use the DOS_initialization function. Will overwrite this function.
             Must be of the type DataInterpolations.LinearInterpolation or a vector of them. It is advised to have extrapolation enabled just in case. 
    - 'egrid': A user made energy grid if not wanting to use build_egrid function. Will overwrite this function. Must be an evenly spaced grid that
               has length(egrid) % 4 == 1 for the numerical integration algorithm to work. 
    - 'BandStructure': Contains splines of [k->E, E->k] for the evaluation of magnetotransport and group velocity

    # Returns
    - The Structure struct with the DOS and egrid assembled or provided by the user
"""
function build_Structure(; las::Laser=build_Laser(), Spatial_DOS::Bool = false, Elemental_System::Int = 1, dimension::Dimension = build_Dimension(),
    bulk_DOS::Union{String,Vector{String},Nothing} = nothing, DOS_folder::Union{String,Vector{String},Nothing} = nothing, 
    bulk_geometry::Union{String,Vector{String},Nothing} = nothing, slab_geometry::Union{String,Vector{String},Nothing} = nothing, 
    atomic_layer_tolerance::Union{Number,Vector{Number}} = 0.1, DOS::Union{spl,Vector{spl},Nothing} = nothing, egrid::Vector{<:Number} = collect(-10.0:0.01:10.0),
    ext_fields = Fields(:(0.0), :(0.0)), bandstructure::Union{spl, Nothing, Vector{spl}} = nothing, FE=0.0)

    DOS = DOS_initialization(bulk_DOS, bulk_geometry, DOS_folder, slab_geometry, atomic_layer_tolerance, dimension, Spatial_DOS, DOS)
    egrid = build_egrid(egrid)
    las_field = calculate_laser_fields(las)
    bandstructure = bandstructure_initialization(bandstructure, DOS, egrid, FE)
    return Structure(Spatial_DOS=Spatial_DOS, Elemental_System=Elemental_System, DOS=DOS, egrid=egrid, dimension=dimension, fields = TotalFields(las_field, ext_fields),
                    bandstructure = bandstructure)
end
"""
    WIP!!!
    DensityMatrix <: SimulationTypes
        Enabled::Bool = false
    end 

    Struct that defines and holds all values for the density matrix propagation.
    This Simulation object doesn't function or couple with the others due to the difference in propagation from coupled 
    ODE to a von-Neumann equation.
"""
@kwdef struct DensityMatrix <: SimulationTypes
    Enabled::Bool
end
"""
    WIP!!!
    build_DensityMatrix(; Enabled = false)

    Once implemented will build a density matrix and store Hamiltonian for propagation via the vonNeumann equation.
"""
function build_DensityMatrix(; Enabled = false)
    return DensityMatrix(Enabled = Enabled)
end
"""
    AthermalElectrons <: SimulationTypes
        Enabled::Bool

        AthermalElectron_ElectronCoupling::Bool # Enables coupling to an electronic bath
        AthermalElectron_PhononCoupling::Bool # Enables coupling to a phononic bath
        Conductivity::Bool # Provides conductivity of a ballistic nature using velocity given by v_g
        EmbeddedAthEM::Bool

        ElectronicRelaxation::Symbol # Implementations are Fermi Liquid Theory (:FLT) or constant (:constant)
        PhononicRelaxation::Symbol # Implementations are constant (:constant) or quasiparticle scattering (:quasi)
        ExcitationMatrixElements::Symbol # Implementation is only match internal energy (:unity)
        Conductive_Velocity::Symbol # Implementation of how gorup velocity is calculated, :constant, :fermigas or :effectiveoneband
        
        FE::Union{Number,Vector{<:Number}} # Shifted Fermi energy to the bottom of the valence band for FLT relaxation and group velocity
        τ::Union{Number,Vector{<:Number}} # Material dependent scale-factor for :FLT relaxation time or the constant value for :constant
        τep::Union{Number,Vector{<:Number}} # Constant relaxation time for phonons
        v_g::Union{Vector{<:Number},Matrix{<:Number}} # Group velocity of electrons calculated assuming a Fermi liquid with μ = FE
    end

    Struct that defines and holds all values for the propagation of athermal electrons
    Enabling this struct assumes an AthEM like system (https://arxiv.org/abs/2503.09479) so can be coupled to electronic
    and phononic thermal baths. Coupling implicitly assumes the other system is enabled. 
"""
@kwdef struct AthermalElectrons <: SimulationTypes
    Enabled::Bool

    AthermalElectron_ElectronCoupling::Bool # Enables coupling to an electronic bath
    AthermalElectron_PhononCoupling::Bool # Enables coupling to a phononic bath
    Conductivity::Bool # Provides conductivity of a ballistic nature using velocity given by v_g
    EmbeddedAthEM::Bool

    ElectronicRelaxation::Symbol # Implementations are Fermi Liquid Theory (:FLT) or constant (:constant)
    PhononicRelaxation::Symbol # Implementations are constant (:constant) or quasiparticle scattering (:quasi)
    ExcitationMatrixElements::Symbol # Implementation is only match internal energy (:unity)
    Conductive_Velocity::Symbol # Implementation of how group velocity is calculated, :constant, :fermigas or :effectiveoneband
    MagnetoTransport::Bool # Whether to add magnetotransport to the problem 
    
    FE::Union{Number,Vector{<:Number}} # Shifted Fermi energy to the bottom of the valence band for FLT relaxation and group velocity
    τ::Union{Number,Vector{<:Number}} # Material dependent scale-factor for :FLT relaxation time or the constant value for :constant
    τep::Union{Number,Vector{<:Number}} # Constant relaxation time for phonons
    v_g::Union{Vector{<:Number},Matrix{<:Number}} # Group velocity of electrons in v(k) calculated assuming a Fermi liquid with μ = FE
end
"""
    build_AthermalElectrons(;structure::Structure, Enabled = false, AthermalElectron_ElectronCoupling = false, 
                            AthermalElectron_PhononCoupling = false, Conductivity = false, ElectronicRelaxation = :FLT, 
                            PhononicRelaxation = :constant, ExcitationMatrixElements = :unity, FE=0.0, τ=1.0, τep = 1000.0, 
                            v_g = nothing, Conductive_Velocity = :constant, EmbeddedAthEM = false)

    Outer constructor function to assemble the AthermalElectrons struct. Unit conversion is detected on all parameters.
    The function will build the group veolcity if one isn't provided by the user.
    Defaults allow any unneccessary parameters for users simulation to be ignored.

    # Arguments
    - 'Enabled': Bool for enabling an athermal electron subssystem
    - 'structure': Structure struct, provide if you want the group velocity calculated for you
    - 'AthermalElectron_ElectronCoupling': Enables athermal electron - thermal electron coupling
    - 'AthermalElectron_PhononCoupling': Enables athermal electron - thermal phonon coupling
    - 'Conductivity': Bool for enabling athermal electron ballistic transport
    - 'ElectronicRelaxation': Method for calculating athermal electron lifetime due to e-e collisions, :FLT or :constant
    - 'PhononicRelaxation': Method for calculating athermal electron lifetime due to e-p collisions, :quasi or :constant
    - 'ExcitationMatrixElements': Decides method to calculate excitation matrix elements, only :unity is currently implemented
    - 'FE': unit = eV: The Fermi energy defined as the difference between the bottom of the valence band in the DOS and 0.0
    - 'τ': unit = fs: A material dependent scalar for the :FLT lifetime or the constant value for :constant e-e lifetime
    - 'τep': unit = fs: The constant lifetime for the athermal electrons due to electron-phonon coupling
    - 'v_g': unit = nm/fs: The group velocity of the ballistic electrons, for the user to define their own group velocity and will overwrite the 
             one calculated by build_group_velocity. Also the value used for a constant velocity.
    - 'Conductive_Velocity': Define the group velocity that build_group_velocity should use, :constant :fermigas, :effectiveoneband
    - 'EmbeddedAthEM': Bool for setting only the surface layer to AthEM and the rest a TTM. Can't be used alongside athermal electron
                       transport

    # Returns
    - The AthermalElectrons struct with the users settings and parameters with any neccessary unit conversion.
"""
function build_AthermalElectrons(; Enabled = false, structure::Structure = build_Structure(), AthermalElectron_ElectronCoupling = false, 
    AthermalElectron_PhononCoupling = false, Conductivity = false, ElectronicRelaxation = :FLT, PhononicRelaxation = :constant, 
    ExcitationMatrixElements = :unity, FE = NaN, τ = NaN, τep = NaN, v_g = nothing, Conductive_Velocity = :constant, EmbeddedAthEM = false,
    MagnetoTransport = false)

    τ = convert_units(u"fs", τ)
    τep = convert_units(u"fs", τep)
    FE = convert_units(u"eV", FE)

    v_g = build_group_velocity(v_g,FE,Conductivity,Conductive_Velocity,structure)

    return AthermalElectrons(Enabled=Enabled, AthermalElectron_ElectronCoupling=AthermalElectron_ElectronCoupling, 
        AthermalElectron_PhononCoupling=AthermalElectron_PhononCoupling, Conductivity=Conductivity, 
        ElectronicRelaxation=ElectronicRelaxation, PhononicRelaxation=PhononicRelaxation, 
        ExcitationMatrixElements=ExcitationMatrixElements, FE=FE, τ=τ, τep=τep, v_g=v_g,
        Conductive_Velocity=Conductive_Velocity,EmbeddedAthEM=EmbeddedAthEM, MagnetoTransport = MagnetoTransport)
end
"""
    ElectronicTemperature <: SimulationTypes
        Enabled::Bool = false

        AthermalElectron_ElectronCoupling::Bool = false # Enables coupling to athermal electrons
        Electron_PhononCoupling::Bool = false # Enables coupling to a phonon thermostat
        Conductivity::Bool = false # Provides diffusive thermal conductivity

        ElectronicHeatCapacity::Symbol = :linear # Whether to use linear (:linear) or non-linear (:nonlinear) 
                                                 # Electronic Heat Capacity
        ElectronPhononCouplingValue::Symbol = :constant # Whether to use constant (:constant) or variable (:variable)
                                                        # electron phonon coupling

        γ::Union{Number,Vector{<:Number}} = 1.0 # Specific heat capacity of electrons at room temperature for linear heat capacity
        κ::Union{Number,Vector{<:Number}} = 1.0 # Thermal conductivity of electrons at room temperature
        λ::Union{Number,Vector{<:Number}} = 1.0 # Electron-phonon mass enhancement factor for non-linear electron-phonon coupling
        ω::Union{Number,Vector{<:Number}} = 1.0 # Second moment of phonon spectral function for non-linear electron-phonon coupling
        g::Union{Number,Vector{<:Number}} = 1.0 # Constant electron-phonon coupling value 
    end

    Struct that defines and holds all values for the propagation of an electronic temperature
    This can be coupled solely to a thermal phonon bath for a Two-Temperature Model simulation or to athermal electrons
    for AthEM propagation with relaxation. Coupling implicitly assumes the other system is enabled.
"""
@kwdef struct ElectronicTemperature <: SimulationTypes
    Enabled::Bool 

    AthermalElectron_ElectronCoupling::Bool # Enables coupling to athermal electrons
    Electron_PhononCoupling::Bool # Enables coupling to a phonon thermostat
    Conductivity::Bool # Provides diffusive thermal conductivity

    ElectronicHeatCapacity::Symbol # Whether to use linear (:linear) or non-linear (:nonlinear) 
                                   # Electronic Heat Capacity
    ElectronPhononCouplingValue::Symbol # Whether to use constant (:constant) or variable (:variable)
                                        # electron phonon coupling

    γ::Union{Number,Vector{<:Number}} # Specific heat capacity of electrons at room temperature for linear heat capacity
    κ::Union{Number,Vector{<:Number}} # Thermal conductivity of electrons at room temperature
    λ::Union{Number,Vector{<:Number}} # Electron-phonon mass enhancement factor for non-linear electron-phonon coupling
    ω::Union{Number,Vector{<:Number}} # Second moment of phonon spectral function for non-linear electron-phonon coupling
    g::Union{Number,Vector{<:Number}} # Constant electron-phonon coupling value 
end
"""
    build_ElectronicTemperature(; Enabled = false, AthermalElectron_ElectronCoupling = false, Electron_PhononCoupling = false, Conductivity = false,
                               ElectronicHeatCapacity = :linear, ElectronPhononCouplingValue = :constant, γ = 1.0, κ = 1.0, λ = 1.0, ω = 1.0, g = 1.0)

    Outer constructor function to assemble the ElectronicTemperature struct. Unit conversion is detected on all parameters.
    Defaults allow any unneccessary parameters for users simulation to be ignored.

    # Arguments
    - 'Enabled': Bool for enabling an thermal electronic bath
    - 'AthermalElectron_ElectronCoupling': Enables athermal electron - thermal electron coupling
    - 'Electron_PhononCoupling': Enables electron bath - phonon bath coupling
    - 'Conductivity': Bool for enabling thermnal electron diffusive transport
    - 'ElectronicHeatCapacity': Method for calculating the electronic heat capacity, :linear or :nonlinear
    - 'ElectronPhononCouplingValue': Method for calculating the electron phonon coupling value, either :constant or :variable
    - 'γ': unit = eV/nm³/K²: Specific heat capacity of electronic bath for :linear ElectronicHeatCapacity
    - 'κ': unit = eV/fs/nm/K: Thermal conductivity of electrons at room temperature
    - 'λ': unit = unitless: Electron-phonon mass enhancement parameter
    - 'ω': unit = eV^2: The second moment of the phonon spectrum
    - 'g': unit = eV/fs/nm³/K: Constant value for the electron-phonon coupling if using :constant

    # Returns
    - The ElectronicTemperature struct with the users settings and parameters with any neccessary unit conversion.
"""
function build_ElectronicTemperature(; Enabled = false, structure=build_Structure(), AthermalElectron_ElectronCoupling = false, Electron_PhononCoupling = false, Conductivity = false,
                               ElectronicHeatCapacity = :linear, ElectronPhononCouplingValue = :constant, γ = NaN, κ = NaN, λ = NaN, ω = NaN, g = NaN)

    γ = convert_units(u"eV/nm^3/K^2", γ)
    if structure.Elemental_System == 1
        κ = convert_units(u"eV/fs/nm/K", κ)
    else
        new_κ = zeros(structure.dimension.length)
        κ = convert_units(u"eV/fs/nm/K", κ)
        for i in eachindex(new_κ)
            X = mat_picker(structure.dimension.grid[i],structure.dimension.InterfaceHeight)
            new_κ[i] = κ[X]
        end
        κ = new_κ
    end
    ω = convert_units(u"eV^2", ω)
    g = convert_units(u"eV/fs/nm^3/K", g)
    return ElectronicTemperature(Enabled=Enabled, AthermalElectron_ElectronCoupling=AthermalElectron_ElectronCoupling, 
                                 Electron_PhononCoupling=Electron_PhononCoupling, Conductivity=Conductivity, 
                                 ElectronicHeatCapacity=ElectronicHeatCapacity, ElectronPhononCouplingValue=ElectronPhononCouplingValue,
                                 γ=γ, κ=κ, λ=λ, ω=ω, g=g)
end
"""
    struct PhononicTemperature <: SimulationTypes
        Enabled::Bool = false

        AthermalElectron_PhononCoupling::Bool = false # Enables coupling to athermal electrons
        Electron_PhononCoupling::Bool = false # Enables coupling to an electron thermostat
        Conductivity::Bool = false # Provides diffusive thermal conductivity

        PhononicHeatCapacity::Symbol = :constant # Whether to use constant (:constant) or non-linear/Simpson's Rule (:nonlinear) 
                                                 # Phononic Heat Capacity
        
        θ::Union{Number,Vector{<:Number}} = 1.0 # Debye temperature for non-linear phonon heat capacity
        n::Union{Number,Vector{<:Number}} = 1.0 # Atomic density for non-linear phonon heat capacity
        Cph::Union{Number,Vector{<:Number}} = 1.0 # Constant phonon heat capacity
        κ::Union{Number,Vector{<:Number}} = 1.0 # Constant phonon thermal conductivity
    end
    Struct that defines and holds all values for the propagation of a phononic temperature
    This can be coupled solely to a thermal electronic bath for a Two-Temperature Model simulation or to athermal electrons
    for AthEM propagation with phonon-relaxation. Coupling implicitly assumes the other system is enabled.
"""
@kwdef struct PhononicTemperature <: SimulationTypes
    Enabled::Bool 

    AthermalElectron_PhononCoupling::Bool # Enables coupling to athermal electrons
    Electron_PhononCoupling::Bool # Enables coupling to an electron thermostat
    Conductivity::Bool # Provides diffusive thermal conductivity

    PhononicHeatCapacity::Symbol # Whether to use constant (:constant) or non-linear/Simpson's Rule (:nonlinear) 
                                 # Phononic Heat Capacity
    
    θ::Union{Number,Vector{<:Number}} # Debye temperature for non-linear phonon heat capacity
    n::Union{Number,Vector{<:Number}} # Atomic density for non-linear phonon heat capacity
    Cph::Union{Number,Vector{<:Number}} # Constant phonon heat capacity
    κ::Union{Number,Vector{<:Number}} # Constant phonon thermal conductivity
end
"""
    build_PhononicTemperature(;Enabled = false, AthermalElectron_PhononCoupling = false, Electron_PhononCoupling = false, 
                               Conductivity = false, PhononicHeatCapacity = :linear, θ = 1.0, n = 1.0, Cph = 1.0, κ = 1.0)

    Outer constructor function to assemble the PhononicTemperature struct. Unit conversion is detected on all parameters.
    Defaults allow any unneccessary parameters for users simulation to be ignored.

    # Arguments
    - 'Enabled': Bool for enabling an thermal electronic bath
    - 'AthermalElectron_PhononCoupling': Enables athermal electron - thermal electron coupling
    - 'Electron_PhononCoupling': Enables electron bath - phonon bath coupling
    - 'Conductivity': Bool for enabling thermal phonon diffusive transport
    - 'PhononicHeatCapacity': Method for calculating the phononic heat capacity, :constant or :nonlinear
    - 'θ': unit = K: Debye temperature of the material
    - 'n': unit = atoms/nm³: Number of atoms per nm³
    - 'Cph': unit = eV/nm³/K: Constant heat capacity for :constant
    - 'κ': unit = eV/nm: Constant thermal conductivity of phonons

    # Returns
    - The PhononicTemperature struct with the users settings and parameters with any neccessary unit conversion.
"""
function build_PhononicTemperature(;Enabled = false, AthermalElectron_PhononCoupling = false, Electron_PhononCoupling = false, 
                                    Conductivity = false, PhononicHeatCapacity = :linear, θ = NaN, n = NaN, Cph = NaN, κ = NaN)

    θ = convert_units(u"K", θ)
    n = convert_units(u"nm^-3", n)
    Cph = convert_units(u"eV/nm^3/K", Cph)
    κ = convert_units(u"eV/nm", κ)
    return PhononicTemperature(Enabled=Enabled, AthermalElectron_PhononCoupling=AthermalElectron_PhononCoupling, 
                               Electron_PhononCoupling=Electron_PhononCoupling, Conductivity=Conductivity, 
                               PhononicHeatCapacity=PhononicHeatCapacity, θ=θ, n=n, Cph=Cph, κ=κ)
end
"""
    W.I.P!!!
    struct ElectronicDistribution <: SimulationTypes
        Enabled::Bool = false

        Electron_PhononCoupling::Bool = false
    end
    Struct that defines and holds all values for the propagation of an electronic distribution
"""
@kwdef struct ElectronicDistribution <: SimulationTypes
    Enabled::Bool = false

    Electron_PhononCoupling::Bool = false
end
"""
    W.I.P!!!
    struct PhononicDistribution <: SimulationTypes
        Enabled::Bool = false

        Electron_PhononCoupling::Bool = false
    end
    Struct that defines and holds all values for the propagation of a phononic distribution
"""
@kwdef struct PhononicDistribution <: SimulationTypes
    Enabled::Bool = false

    Electron_PhononCoupling::Bool = false
end
"""
    struct Simulation <: SimulationTypes
        densitymatrix::DensityMatrix
        electronictemperature::ElectronicTemperature
        phononictemperature::PhononicTemperature
        athermalelectrons::AthermalElectrons
        electronicdistribution::ElectronicDistribution
        phononicdistribution::PhononicDistribution
        structure::Structure
        laser::Laser
    end
    This struct contains all the others and is the main simulation object both in assembling a simulation and during it
"""
@kwdef struct Simulation <: SimulationTypes
    densitymatrix::DensityMatrix
    electronictemperature::ElectronicTemperature
    phononictemperature::PhononicTemperature
    athermalelectrons::AthermalElectrons
    electronicdistribution::ElectronicDistribution
    phononicdistribution::PhononicDistribution
    structure::Structure
    laser::Laser
end
"""
    build_Simulation(;densitymatrix::Union{DensityMatrix,NamedTuple,Nothing}=nothing, electronictemperature::Union{ElectronicTemperature,NamedTuple,Nothing}=nothing,
                           phononictemperature::Union{PhononicTemperature,NamedTuple,Nothing}=nothing, athermalelectrons::Union{AthermalElectrons,NamedTuple,Nothing}=nothing,
                           electronicdistribution::Union{ElectronicDistribution,NamedTuple,Nothing}=nothing, phononicdistribution::Union{PhononicDistribution,NamedTuple,Nothing}=nothing,
                           structure::Union{Structure,NamedTuple,Nothing}=nothing, laser::Union{Laser,NamedTuple,Nothing}=nothing)

    Assembles the full Simulation struct from the requested components. Any systems not provided to the function are disabled by default.
    The user can send either a completed Struct of the correct type or a NamedTuple with the correct key-word arguments to assemble the struct directly within the function
    using the representative build_x function where x is the subsystem.

    # Arguments
    - 'densitymatrix': The DensityMatrix subsystem
    - 'electronictemperature': The ElectronicTemperature subsystem
    - 'phononictemperature': The PhononicTemperature subsystem
    - 'athermalelectrons': The AthermalElectrons subsystem
    - 'electronicdistribution': The ElectronicDistribution subsystem
    - 'phononicdistribution': The PhononicDistribution subsystem
    - 'structure': The Structure subsystem
    - 'laser': The Laser subsystem

    # Returns
    - The PhononicTemperature struct with the users settings and parameters with any neccessary unit conversion.
"""
function build_Simulation(;densitymatrix::Union{DensityMatrix,NamedTuple,Nothing}=nothing, electronictemperature::Union{ElectronicTemperature,NamedTuple,Nothing}=nothing,
                           phononictemperature::Union{PhononicTemperature,NamedTuple,Nothing}=nothing, athermalelectrons::Union{AthermalElectrons,NamedTuple,Nothing}=nothing,
                           electronicdistribution::Union{ElectronicDistribution,NamedTuple,Nothing}=nothing, phononicdistribution::Union{PhononicDistribution,NamedTuple,Nothing}=nothing,
                           structure::Union{Structure,NamedTuple,Nothing}=nothing, laser::Union{Laser,NamedTuple,Nothing}=nothing)

    #Implement calculation of energy grid and v_g for athermal transport
    temp = NamedTuple()

    densitymatrix = if densitymatrix isa DensityMatrix
        densitymatrix
    else
        build_DensityMatrix(; merge(temp, densitymatrix isa NamedTuple ? densitymatrix : NamedTuple())...)
    end

    electronictemperature = if electronictemperature isa ElectronicTemperature
        electronictemperature
    else
        build_ElectronicTemperature(; merge(temp, electronictemperature isa NamedTuple ? electronictemperature : NamedTuple())...)
    end

    phononictemperature = if phononictemperature isa PhononicTemperature
        phononictemperature
    else
        build_PhononicTemperature(; merge(temp, phononictemperature isa NamedTuple ? phononictemperature : NamedTuple())...)
    end

    athermalelectrons = if athermalelectrons isa AthermalElectrons
        athermalelectrons
    else
        build_AthermalElectrons(; merge(temp, athermalelectrons isa NamedTuple ? athermalelectrons : NamedTuple())...)
    end

    electronicdistribution = if electronicdistribution isa ElectronicDistribution
        electronicdistribution
    else
        ElectronicDistribution(; merge(temp, electronicdistribution isa NamedTuple ? electronicdistribution : NamedTuple())...)
    end

    phononicdistribution = if phononicdistribution isa PhononicDistribution
        phononicdistribution
    else
        PhononicDistribution(; merge(temp, phononicdistribution isa NamedTuple ? phononicdistribution : NamedTuple())...)
    end

    structure = if structure isa Structure
        structure
    else
        build_Structure(; merge(temp, structure isa NamedTuple ? structure : NamedTuple())...)
    end

    laser = if laser isa Laser
        laser
    else
        build_Laser(; merge(temp, laser isa NamedTuple ? laser : NamedTuple())...)
    end

    return Simulation(densitymatrix=densitymatrix,electronictemperature=electronictemperature,phononictemperature=phononictemperature,athermalelectrons=athermalelectrons,
    electronicdistribution=electronicdistribution,phononicdistribution=phononicdistribution,structure=structure,laser=laser)
end
"""
    build_egrid(hv)

    Builds an energy grid with spacing of 0.01 eV, limits of -2*hv to 2*hv and ensures that length(egrid) % 4 == 1

    # Arguments
    - 'hv': Photon energy of the laser

    # Returns
    - Evenly spaced energy grid that is suitable for the numerical integration algorithm and of sufficient accuracy/discretisation for accurate dynamics
"""
function build_egrid(egrid)
    side = 1
    dE = egrid[2] - egrid[1]
    while length(egrid) % 4 != 1
        if side == 1
            push!(egrid,egrid[end]+dE)
            side = -1
        elseif side == -1
            pushfirst!(egrid,egrid[1]-dE)
            side = 1
        end
    end
    return egrid
end
"""
    SimulationTypes
    Abstract type that holds all information for a SimulationTypes
"""
abstract type SimulationTypes end
"""
    spl=DataInterpolations.LinearInterpolation
    A convenience type definition to make type specificity easier throughout the code
"""
global const spl=DataInterpolations.LinearInterpolation
"""
    Constants = (ħ = 8.617e-5,kB = 0.6582,me = 3.109e-31)
    Global named tuple for accessing constant physical values during a Simulation
"""
global const Constants = (ħ = convert_units(Unitful.ħ),kB = convert_units(Unitful.k),me = convert_units(Unitful.me))
"""
    WIP!!!
    struct DensityMatrix <: SimulationTypes
        Enabled::Bool = false
    end 
    Struct that defines and holds all values for the density matrix propagation.
    This Simulation object doesn't function or couple with the others due to the difference in propagation from coupled 
    ODE to a von-Neumann equation.
"""
@kwdef struct DensityMatrix <: SimulationTypes
    Enabled::Bool
end

function DensityMatrix(;Enabled=false)
    return DensityMatrix(Enabled=Enabled)
end
"""
    struct AthermalElectrons <: SimulationTypes
        Enabled::Bool = false

        AthermalElectron_ElectronCoupling::Bool = false #Enables coupling to an electronic bath
        AthermalElectron_PhononCoupling::Bool = false #Enables coupling to a phononic bath
        Conductivity::Bool = false #Provides conductivity of a ballistic nature using velocity given by v_g

        ElectronicRelaxation::Symbol = :FLT #Implementations are Fermi Liquid Theory (:FLT) or constant (:constant)
        PhononicRelaxation::Symbol = :constant #Implementations are constant (:constant) or quasiparticle scattering (:quasi)
        ExcitationMatrixElements::Symbol = :laser_equal_internal_E #Implementation is only match internal energy (:laser_equal_internal_E)
        
        FE::Union{Real,Vector{<:Real}} = 0.0 #Shifted Fermi energy to the bottom of the valence band for FLT relaxation and group velocity
        τ::Union{Real,Vector{<:Real}} = 1.0 #Material dependent scale-factor for FLT relaxation time
        τee::Union{Real,Vector{<:Real}} = 1.0 #Constant relaxation time for electrons
        τep::Union{Real,Vector{<:Real}} = 1000.0 #Constant relaxation time for phonons
        v_g::Union{Vector{Real},Vector{Vector{Real}}} = [1.0] #Group velocity of electrons calculated assuming a Fermi liquid with μ = FE
    end
    Struct that defines and holds all values for the propagation of athermal electrons
    Enabling this struct assumes an AthEM like system (https://arxiv.org/abs/2503.09479) so can be coupled to electronic
    and phononic thermal baths. Coupling implicitly assumes the other system is enabled. 
"""
@kwdef struct AthermalElectrons <: SimulationTypes
    Enabled::Bool

    AthermalElectron_ElectronCoupling::Bool #Enables coupling to an electronic bath
    AthermalElectron_PhononCoupling::Bool #Enables coupling to a phononic bath
    Conductivity::Bool #Provides conductivity of a ballistic nature using velocity given by v_g
    EmbeddedAthEM::Bool

    ElectronicRelaxation::Symbol #Implementations are Fermi Liquid Theory (:FLT) or constant (:constant)
    PhononicRelaxation::Symbol #Implementations are constant (:constant) or quasiparticle scattering (:quasi)
    ExcitationMatrixElements::Symbol #Implementation is only match internal energy (:unity)
    Conductive_velocity::Symbol #Implementation of how gorup velocity is calculated, :constant, :fermigas or :effectiveoneband
    
    FE::Union{Real,Vector{<:Real}} #Shifted Fermi energy to the bottom of the valence band for FLT relaxation and group velocity
    τ::Union{Real,Vector{<:Real}} #Material dependent scale-factor for :FLT relaxation time or the constant value for :constant
    τee::Union{Real,Vector{<:Real}} #Constant relaxation time for electrons
    τep::Union{Real,Vector{<:Real}} #Constant relaxation time for phonons
    v_g::Union{Vector{<:Real},Vector{Vector{<:Real}}} #Group velocity of electrons calculated assuming a Fermi liquid with μ = FE
end

function AthermalElectrons(;Enabled=false, AthermalElectron_ElectronCoupling=false, AthermalElectron_PhononCoupling=false, Conductivity=false,
                           ElectronicRelaxation=:FLT, PhononicRelaxation=:constant, ExcitationMatrixElements=:unity,
                           FE=0.0, τ=1.0, τee=1.0, τep=1000.0, v_g=[1.0], conductive_velocity=:constant)
    
    FE = convert_units(FE)
    τ = convert_units(τ)
    τee = convert_units(τee)
    τep = convert_units(τep)
    v_g = convert_units(v_g)
    return AthermalElectrons(Enabled=Enabled, AthermalElectron_ElectronCoupling=AthermalElectron_ElectronCoupling, 
                             AthermalElectron_PhononCoupling=AthermalElectron_PhononCoupling, Conductivity=Conductivity, 
                             ElectronicRelaxation=ElectronicRelaxation, PhononicRelaxation=PhononicRelaxation, 
                             ExcitationMatrixElements=ExcitationMatrixElements, Fe=Fe, τ=τ, τee=τee, τep=τep, v_g=v_g,
                             conductive_velocity=conductive_velocity)
end
"""
    struct ElectronicTemperature <: SimulationTypes
        Enabled::Bool = false

        AthermalElectron_ElectronCoupling::Bool = false #Enables coupling to athermal electrons
        Electron_PhononCoupling::Bool = false #Enables coupling to a phonon thermostat
        Conductivity::Bool = false #Provides diffusive thermal conductivity

        ElectronicHeatCapacity::Symbol = :linear #Whether to use linear (:linear) or non-linear (:nonlinear) 
                                            #Electronic Heat Capacity
        ElectronPhononCouplingValue::Symbol = :constant #Whether to use constant (:constant) or variable (:variable)
                                                        #electron phonon coupling

        γ::Union{Real,Vector{<:Real}} = 1.0 #Specific heat capacity of electrons at room temperature for linear heat capacity
        κ::Union{Real,Vector{<:Real}} = 1.0 #Thermal conductivity of electrons at room temperature
        λ::Union{Real,Vector{<:Real}} = 1.0 #Electron-phonon mass enhancement factor for non-linear electron-phonon coupling
        ω::Union{Real,Vector{<:Real}} = 1.0 #Second moment of phonon spectral function for non-linear electron-phonon coupling
        g::Union{Real,Vector{<:Real}} = 1.0 #Constant electron-phonon coupling value 
    end
    Struct that defines and holds all values for the propagation of an electronic temperature
    This can be coupled solely to a thermal phonon bath for a Two-Temperature Model simulation or to athermal electrons
    for AthEM propagation with relaxation. Coupling implicitly assumes the other system is enabled.
"""
@kwdef struct ElectronicTemperature <: SimulationTypes
    Enabled::Bool 

    AthermalElectron_ElectronCoupling::Bool #Enables coupling to athermal electrons
    Electron_PhononCoupling::Bool #Enables coupling to a phonon thermostat
    Conductivity::Bool #Provides diffusive thermal conductivity

    ElectronicHeatCapacity::Symbol #Whether to use linear (:linear) or non-linear (:nonlinear) 
                                   #Electronic Heat Capacity
    ElectronPhononCouplingValue::Symbol #Whether to use constant (:constant) or variable (:variable)
                                        #electron phonon coupling

    γ::Union{Real,Vector{<:Real}} #Specific heat capacity of electrons at room temperature for linear heat capacity
    κ::Union{Real,Vector{<:Real}} #Thermal conductivity of electrons at room temperature
    λ::Union{Real,Vector{<:Real}} #Electron-phonon mass enhancement factor for non-linear electron-phonon coupling
    ω::Union{Real,Vector{<:Real}} #Second moment of phonon spectral function for non-linear electron-phonon coupling
    g::Union{Real,Vector{<:Real}} #Constant electron-phonon coupling value 
end

function ElectronicTemperature(;Enabled=false, AthermalElectron_ElectronCoupling=false, Electron_PhononCoupling=false, Conductivity=false,
                               ElectronicHeatCapacity = :linear, ElectronPhononCouplingValue=:constant, γ=1.0, κ=1.0, λ=1.0, ω=1.0, g=1.0)

    γ = convert_units(γ)
    κ = convert_units(κ)
    λ = convert_units(λ)
    ω = convert_units(ω)
    g = convert_units(g)
    return ElectronicTemperature(Enabled=Enabled, AthermalElectron_ElectronCoupling=AthermalElectron_ElectronCoupling, 
                                 Electron_PhononCoupling=Electron_PhononCoupling, Conductivity=Conductivity, 
                                 ElectronicHeatCapactity=ElectronicHeatCapacity, ElectronPhononCouplingValue=ElectronPhononCouplingValue,
                                 γ=γ, κ=κ, λ=λ, ω=ω, g=g)
end
"""
    struct PhononicTemperature <: SimulationTypes
        Enabled::Bool = false

        AthermalElectron_PhononCoupling::Bool = false #Enables coupling to athermal electrons
        Electron_PhononCoupling::Bool = false #Enables coupling to an electron thermostat
        Conductivity::Bool = false #Provides diffusive thermal conductivity

        PhononicHeatCapacity::Symbol = :constant #Whether to use constant (:constant) or non-linear/Simpson's Rule (:nonlinear) 
                                                #Phononic Heat Capacity
        
        θ::Union{Real,Vector{<:Real}} = 1.0 #Debye temperature for non-linear phonon heat capacity
        n::Union{Real,Vector{<:Real}} = 1.0 #Atomic density for non-linear phonon heat capacity
        Cph::Union{Real,Vector{<:Real}} = 1.0 #Constant phonon heat capacity
        κ::Union{Real,Vector{<:Real}} = 1.0 #Constant phonon thermal conductivity
    end
    Struct that defines and holds all values for the propagation of a phononic temperature
    This can be coupled solely to a thermal electronic bath for a Two-Temperature Model simulation or to athermal electrons
    for AthEM propagation with phonon-relaxation. Coupling implicitly assumes the other system is enabled.
"""
@kwdef struct PhononicTemperature <: SimulationTypes
    Enabled::Bool 

    AthermalElectron_PhononCoupling::Bool #Enables coupling to athermal electrons
    Electron_PhononCoupling::Bool #Enables coupling to an electron thermostat
    Conductivity::Bool #Provides diffusive thermal conductivity

    PhononicHeatCapacity::Symbol #Whether to use constant (:constant) or non-linear/Simpson's Rule (:nonlinear) 
                                 #Phononic Heat Capacity
    
    θ::Union{Real,Vector{<:Real}} #Debye temperature for non-linear phonon heat capacity
    n::Union{Real,Vector{<:Real}} #Atomic density for non-linear phonon heat capacity
    Cph::Union{Real,Vector{<:Real}} #Constant phonon heat capacity
    κ::Union{Real,Vector{<:Real}} #Constant phonon thermal conductivity
end

function PhononicTemperature(;Enabled=false, AthermalElectron_PhononCoupling=false, Electron_PhononCoupling=false, Conductivity=false,
                             PhononicHeatCapacity = :linear, θ=1.0, n=1.0, Cph=1.0, κ=1.0)

    θ = convert_units(θ)
    n = convert_units(n)
    Cph = convert_units(Cph)
    κ = convert_units(κ)
    return PhononicTemperature(Enabled=Enabled, AthermalElectron_PhononCoupling=AthermalElectron_PhononCoupling, 
                               Electron_PhononCoupling=Electron_PhononCoupling, Conductivity=Conductivity, 
                               PhononicHeatCapactity=PhononicHeatCapacity, θ=θ, n=n, Cph=Cph, κ=κ)
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
    struct ElectronicDistribution <: SimulationTypes
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
    struct Dimension <: SimulationTypes
        length::Union{Int,Vector{Int}}=1 #The length of the grid, not the depth of the slab
        grid::AbsractArray{<:Real}=[0.0] #The grid the simulation is solved over
        spacing::Union{Real,Vector<:Real}=0.0 #The spacing between grid points
    end
    Struct that contains all information regarding the grid that the simulation is performed on.
"""
@kwdef struct Dimension <: SimulationTypes
    length::Union{Int,Vector{Int}} #The length of the grid, not the depth of the slab
    grid::AbstractArray{<:Real} #The grid the simulation is solved over
    spacing::Union{Real,Vector{<:Real}} #The spacing between grid points
end

function Dimension(grid=[0.0])
    L = length(grid)
    grid = convert_units(grid)
    if L > 1
        spacing = grid[2]-grid[1]
    else
        spacing = 1.0
    end
    return Dimension(length=L, grid=grid, spacing=spacing)
end
"""
    struct Structure <: SimulationTypes
        Spatial_DOS::Bool = false #Whether to vary the DOS with height - if so the DOS becomes a vector

        Elemental_System::Int = 1 #The number of elemental systems, if < 1 then each constant and vector
                                  #of material parameters needs to become a vector of length=Elemental_System

        DOS::Union{spl,Vector{spl}} = LinearInterpolation([5,6,7,8],[1,2,3,4]) #The density of states of the simulation
        egrid::Vector{Real} = [0.0] #An energy grid for electronic or phononic distributions to be solved on

        dimension::Dimension = Homogenous() #A struct holding all spatial grid structure (0D or 1D)
    end
    Struct that contains any spatial information including the DOS, the spatial grid to solve the problem on and
    the elemental composition of the simulation (e.g. an antenna-reactor system would contain two elemental systems)
"""
@kwdef struct Structure <: SimulationTypes
    Spatial_DOS::Bool #Whether to vary the DOS with height - if so the DOS becomes a vector

    Elemental_System::Int #The number of elemental systems, if < 1 then each constant and vector
                          #of material parameters needs to become a vector of length=Elemental_System

    DOS::Union{spl,Vector{spl}} #The density of states of the simulation
    egrid::Vector{Real} #An energy grid for electronic or phononic distributions to be solved on

    dimension::Dimension #A struct holding all spatial grid structure (0D or 1D)
end

function Structure(;Spatial_DOS=false, Elemental_System=1, DOS=LinearInterpolation([5,6,7,8],[1,2,3,4]), egrid=[0.0], 
                   dimension=Dimension())
    
    egrid = convert_units(egrid)

    return Structure(Spatial_DOS=Spatial_DOS, Elemental_System=Elemental_System, DOS=DOS, egrid=egrid, dimension=dimension)
end
"""
    struct Laser <: SimulationTypes
        envelope::Symbol = :Gaussian #Currently implemented are :Gaussian, :HyperbolicSecant, :Lorentzian and :Rectangular
        FWHM::Real = 10.0 #The Full-Width Half-Maximum of the laser, for rectnagular half the length
        Power::Real = 10.0 #The unabsorbed fluence of the laser
        hv::Real = 1.55 #The photon frequency of the laser
        Transport::Symbol = :Optical #:Optical, :Ballistic and :Combined are the options for how the laser decays into a slab
        ϵ::Union{Real,Vector{<:Real}} = 1.0 #The inverse of the absorption coefficient
        R::Real = 0.0 #The reflectivity of the sample
        δb::Union{Real,Vector{<:Real}} = 1.0 #The ballistic length of electrons
    end
    Struct that contains the shape and parameters that define the laser, these also include some material parameters
    such as the inverse absorption coefficient (ϵ) and reflectivity (R).
"""
@kwdef struct Laser <: SimulationTypes
    envelope::Symbol = :Gaussian #Currently implemented are :Gaussian, :HyperbolicSecant, :Lorentzian and :Rectangular
    Transport::Symbol #:Optical, :Ballistic and :Combined are the options for how the laser decays into a slab

    FWHM::Real #The Full-Width Half-Maximum of the laser, for rectnagular half the length
    ϕ::Real #The unabsorbed fluence of the laser
    hv::Real #The photon frequency of the laser
    ϵ::Union{Real,Vector{<:Real}} #The inverse of the absorption coefficient
    R::Real #The reflectivity of the sample
    δb::Union{Real,Vector{<:Real}} #The ballistic length of electrons
end

function Laser(;envelope=:Gaussian, FWHM=10.0, Power=10.0, hv=5.0, Transport=:Optical, ϵ=1.0, R=0.0, δb=1.0)
    
    FWHM = convert_units(FWHM)
    Power = convert_units(Power)
    hv = convert_units(hv)
    ϵ = convert_units(ϵ)
    δb = convert_units(δb)
    return Laser(envelope=envelope, FWHM=FWHM, Power=Power, hv=hv, Transport=Transport, ϵ=ϵ, R=R, δb=δb)
end
"""
    struct Systems <: SimulationTypes
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

function build_Simulation(;densitymatrix::Union{DensityMatrix,NamedTuple,Nothing}=nothing,electronictemperature::Union{ElectronicTemperature,NamedTuple,Nothing}=nothing,
                           phononictemperature::Union{PhononicTemperature,NamedTuple,Nothing}=nothing,athermalelectrons::Union{AthermalElectrons,NamedTuple,Nothing}=nothing,
                           electronicdistribution::Union{ElectronicDistribution,NamedTuple,Nothing}=nothing,phononicdistribution::Union{PhononicDistribution,NamedTuple,Nothing}=nothing,
                           structure::Union{Structure,NamedTuple,Nothing}=nothing,laser::Union{Laser,NamedTuple,Nothing}=nothing,
                           bulk_DOS::Union{String,Vector{String},Nothing}=nothing,DOS_folder::Union{String,Vector{String},Nothing}=nothing,
                           bulk_geometry::Union{String,Vector{String},Nothing}=nothing,slab_geometry::Union{String,Vector{String},Nothing}=nothing,atomic_layer_tolerance::Union{Real,Vector{Real}}=0.1)

    #Implement calculation of energy grid and v_g for athermal transport
    structure_tuple,athermal_tuple = DOS_FE_initialization(bulk_DOS,bulk_geometry,DOS_folder,slab_geometry,structure,atomic_layer_tolerance)
    grid = build_energy_grid(laser,structure)
    group_velocity = build_group_velocity(athermalelectrons,structure,grid,athermal_tuple,dos_tuple)

    structure_tuple=merge(structure_tuple,grid)
    athermal_tuple=merge(athermal_tuple,group_velocity)

    temp = NamedTuple()

    densitymatrix = if densitymatrix isa DensityMatrix
        densitymatrixS
    else
        DensityMatrix(; merge(temp, densitymatrix isa NamedTuple ? densitymatrix : NamedTuple())...)
    end

    electronictemperature = if electronictemperature isa ElectronicTemperature
        electronictemperature
    else
        ElectronicTemperature(; merge(temp, electronictemperature isa NamedTuple ? electronictemperature : NamedTuple())...)
    end

    phononictemperature = if phononictemperature isa PhononicTemperature
        phononictemperature
    else
        PhononicTemperature(; merge(temp, phononictemperature isa NamedTuple ? phononictemperature : NamedTuple())...)
    end

    athermalelectrons = if athermalelectrons isa AthermalElectrons
        athermalelectrons
    else
        AthermalElectrons(; merge(temp,athermal_tuple, athermalelectrons isa NamedTuple ? athermalelectrons : NamedTuple())...)
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
        Structure(; merge(temp,structure_tuple, structure isa NamedTuple ? structure : NamedTuple())...)
    end

    laser = if laser isa Laser
        laser
    else
        Laser(; merge(temp, laser isa NamedTuple ? laser : NamedTuple())...)
    end

    return Systems(densitymatrix=densitymatrix,electronictemperature=electronictemperature,phononictemperature=phononictemperature,athermalelectrons=athermalelectrons,
    electronicdistribution=electronicdistribution,phononicdistribution=phononicdistribution,structure=structure,laser=laser)
end
"""
    There is no more dim so need to organise that as well as needs to come from structure.dimension.length as well as number of atoms 
"""
function DOS_FE_initialization(bulk_DOS,bulk_geometry,DOS_folder,slab_geometry,structure,atomic_layer_tolerance)
    length = isnothing(structure) ? 1 : structure.length
    if !isnothing(bulk_DOS)
        if length(bulk_DOS) == 1
            FE=get_FermiEnergy(bulk_DOS)
            Vbulk = get_unitcellvolume(bulk_geometry)
            if structure.spacing == true 
                dim = isnothing(structure) ? Dimension() : structure.dimension
                DOS = spatial_DOS(DOS_folder,slab_geometry,bulk_DOS,Vbulk,dim,layer_tolerance)
            else
                DOS = generate_DOS(bulk_DOS,1/Vbulk) 
            end
            dos_tuple = (DOS=DOS,)
            athermal_tuple =(FE=FE,)
        else
            FE_vec = zeros(length(bulk_DOS))
            Vbulk_vec = zeros(length(bulk_DOS))
            DOS_vec = Vector(undef,length(bulk_DOS))
            for i in eachindex(bulk_DOS)
                FE_vec[i] = get_FermiEnergy(bulk_DOS[i])
                Vbulk_vec[i] = get_unitcellvolume(bulk_geometry[i])
                if structure.spacing == true 
                    dim = isnothing(structure) ? Dimension() : structure.dimension
                    DOS_vec[i] = spatial_DOS(DOS_folder[i],slab_geometry[i],bulk_DOS[i],Vbulk[i],dim,atomic_layer_tolerance)
                else
                    DOS_vec[i] = generate_DOS(bulk_DOS[i],1/Vbulk[i])
                end
            end
            dos_tuple = (DOS=DOS_vec,)
            athermal_tuple =(FE=FE_vec,)
        end
    else
        dos_tuple = NamedTuple()
        athermal_tuple = NamedTuple()
    end
    return dos_tuple,athermal_tuple
end

function build_energy_grid(laser,structure)
    if structure != nothing && :egrid in propertynames(structure)
        egrid = structure.egrid
        side = 1
        while length(egrid) % 4 != 1
            if side == 1
                push!(structure.egrid,egrid[end]+(egrid[end]-egrid[end-1]))
                side = -1
            elseif side == -1
                pushfirst!(structure.egrid,egrid[1]-(egrid[2]-egrid[1]))
                side = 1
            end
        end
        return NamedTuple()
    else
        if :hv in propertynames(laser)
            hv = laser.hv
            return (egrid=build_egrid(hv),)
        else
            hv = Laser().hv
            return (egrid=build_egrid(hv),)
        end

    end
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

function build_group_velocity(athermalelectrons,structure,grid,athermal_tuple,dos_tuple)
    if :conductivity in propertynames(athermalelectrons) && athermalelectrons.conductivity == true
        if athermalelectrons.Conductive_velocity == :constant
            return NamedTuple()
        elseif athermalelectrons.Conductive_velocity == :fermigas
            egrid = if grid == NamedTuple()
                structure.egrid
            else
                grid.egrid
            end
            EF = if athermal_tuple == NamedTuple()
                athermalelectrons.FE
            else
                athermal_tuple.FE
            end
            return (v_g = get_fermigas_velocity.(Ref(egrid),EF),)
        elseif athermalelectrons.Conductive_velocity == :effectiveoneband
            egrid = if grid == NamedTuple()
                structure.egrid
            else
                grid.egrid
            end
            EF = if athermal_tuple == NamedTuple()
                athermalelectrons.FE
            else
                athermal_tuple.FE
            end
            DOS = if dos_tuple == NamedTuple()
                structure.DOS
            else
                dos_tuple.DOS
            end

            return (v_g = effective_one_band_velocity(DOS,egrid,EF),)
        end
    else
        return NamedTuple()
    end 
end

function get_fermigas_velocity(egrid,EF)
    return sqrt.(2*(egrid.+EF)./Constants.me)
end

function effective_one_band_velocity(DOS,egrid,FE)
    k_E = get_dispersionrelation(DOS,egrid,FE)
    v_g = similar(k_E)
    if eltype(v_g) <: AbstractVector
        for j in eachindex(v_g)
            for i in eachindex(k_E[j])
                ΔE = egrid[2]-egrid[1]
                if i == 1 
                    v_g[j][i] = 1/Constants.ħ * 1/((k_E[j][i+1] - k_E[j][i]) / ΔE)
                elseif i == length(k_E)
                    v_g[j][i] = 1/Constants.ħ * 1/((k_E[j][i] - k_E[j][i-1]) / ΔE)
                else
                    cd = (k_E[i+1]-k_E[i-1]) / (2*ΔE)
                    v_g[j][i] = 1/Constants.ħ * (1 / cd)
                end
            end
        end
    else
        for i in eachindex(k_E)
            ΔE = egrid[2]-egrid[1]
            if i == 1 
                v_g[i] = 1/Constants.ħ * 1/((k_E[i+1] - k_E[i]) / ΔE)
            elseif i == length(k_E)
                v_g[i] = 1/Constants.ħ * 1/((k_E[i] - k_E[i-1]) / ΔE)
            else
                cd = (k_E[i+1]-k_E[i-1]) / (2*ΔE)
                v_g[i] = 1/Constants.ħ * (1 / cd)
            end
        end
    end
    return v_g
end

function get_dispersionrelation(DOS,egrid,FE)
    if DOS isa Vector{spl}
        k_E=fill(zeros(egrid),length(DOS))
    else
        k_E = zeros(length(egrid))
    end
    
    factor = 6*pi #6π^2/σ where σ is a spin factor (2 for electrons)
    int(u,p) = DOS(u)

    if DOS isa Vector{spl}
        for v in eachindex(k_E)
            for (i,E) in enumerate(egrid)
                prob=IntegralProblem(int,-FE,E)
                sol = solve(prob,HCubatureJL(initdiv=100),abstol=1e-8,reltol=1e-8)
                k_E[v][i] = cbrt(factor*sol.u)
            end
        end
    else 
        for (i,E) in enumerate(egrid)
            prob=IntegralProblem(int,-FE,E)
            sol = solve(prob,HCubatureJL(initdiv=100),abstol=1e-8,reltol=1e-8)
            k_E[i] = cbrt(factor*sol.u)
        end
    end
    return k_E
end
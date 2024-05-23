@kwdef struct MaterialParameters <: simulation
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
    DOSFile::String #File location of the DOS data
    λ::Float64 #Second momentum of spectral function
    g::Float64 # Linear Electron-phonon coupling constant
    ballistic::Real # Ballistic length of electrons
end

function define_material_parameters(;extcof,fermien,gamma,debye,noatoms,
    plasma,thermalcond,elecperatom,eleceffmass,dos,secmomspecfun,elecphon)

    matpat=MaterialParameters(ϵ=extcof,FE=fermien,γ=gamma
    ,θ=debye,n=noatoms,ω=plasma,κ=thermalcond,ne=elecperatom,
    effmass=eleceffmass,DOSFile=dos,λ=secmomspecfun,g=elecphon)

    return matpat
end

function define_material_parameters(dict::Dict;extcof=dict["ExtCof"],fermien=dict["FE"]
    ,gamma=dict["Gamma"],debye=dict["Debye"],noatoms=dict["AtomDens"],plasma=dict["Plasma"]
    ,thermalcond=dict["RTKappa"],elecperatom=dict["ne"],eleceffmass=dict["EffMass"]
    ,dos=dict["DOS"],secmomspecfun=Dict["SpectralFunc"],elecphon=Dict["g"])

    matpat=MaterialParameters(ϵ=extcof,FE=fermien,γ=gamma
    ,θ=debye,n=noatoms,ω=plasma,κ=thermalcond,ne=elecperatom,effmass=eleceffmass,
    DOSFile=dos,λ=secmomspecfun,g=elecphon)

    return matpat
end
module Structs

using Unitful,DelimitedFiles,StaticArrays

@kwdef struct GlobalConsts
    kB::Float64 #Boltzmann Constant
    elecmass::Float64 #Electron mass
    eleccharge::Float64 #Electronic charge
    ϵ0::Float64 #Dielectric constant in vacuum
    hbar::Float64 #Reduced plank constant
end

@kwdef struct LaserParam
    FWHM::Float64 #Full-Width at Half-Maximum
    ϕ::Float64 #Fluence
    delay::Float64 #Peak offset from 0
    hv::Float64 #Photon energy
end

@kwdef struct MaterialParam
    α::Float64 #Absorption Coefficient
    FE::Float64 #Fermi Energy
    γ::Float64 #Specific heat capacity
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
end

@kwdef mutable struct SimMaterialParam
    Tel::Vector{Float64} #Electronic temperature
    Tph::Vector{Float64} #Phononic temperature
    μ::Float64 #Chemical Potential
    zgrid::Vector{Float64}
    DOS #Interpolate Density of States
    lb::Float64 #lower bound of the DOS data
    ub::Float64 #Upper bound of the DOS data
    NumElec::Float64 #Number of electrons defined by integral of f(E,μ,T)*D(E)
    τ::Float64 #Material scaling of the e^*-e relaxation time
    g::Float64 #Electron-phonon coupling constant
end

@kwdef struct SimParam
    ChemPot::String
    ElecPhon::String
    ElecHeat::String
    Dim::Int16
    Output::String
    vers::String
    SimEnd::Float64
end

function Units(Input)
    Input["g"]=ustrip(uconvert(u"eV/nm^3/fs/K", Input["g"]u"W/m^3/K"))
    Input["gamma"]=ustrip(uconvert(u"eV/nm^3/K^2",Input["gamma"]u"J/m^3/K^2"))
    Input["ExtCof"]=ustrip(uconvert(u"nm",Input["ExtCof"]u"m"))
    Input["AtomDens"]=ustrip(uconvert(u"nm^-3",Input["AtomDens"]u"m^-3"))
    Input["FWHM"]=ustrip(uconvert(u"fs",Input["FWHM"]u"s"))
    Input["Fluence"]=ustrip(uconvert(u"eV/nm^2",Input["Fluence"]u"J/m^2"))
    Input["LaserOff"]=ustrip(uconvert(u"fs",Input["LaserOff"]u"s"))
    Input["RTKappa"]=ustrip(uconvert(u"eV/nm/K/fs",Input["RTKappa"]u"W/m/K"))
    Input["SimEnd"]=ustrip(uconvert(u"fs",Input["SimEnd"]u"s"))
    Input["Length"]=ustrip(uconvert(u"nm",Input["Length"]u"m"))
    Input["dz"]=ustrip(uconvert(u"nm",Input["dz"]u"m"))
    Input["Plasma"]=Input["Plasma"]/(6.582e-1*2*pi)
    return Input
end

function parameterbuilder(InputFile::String)
    IO = readdlm(InputFile,'=',comments=true,comment_char='#')
    for (x,y) in enumerate(IO[:,2])
        if typeof(y) == SubString{String}
            IO[x,2]=strip(IO[x,2],' ')
        end
    end
    Input=Dict(IO[i,1]=>IO[i,2] for i in 1:size(IO,1))
    if Input["Type"]=="E2TM"
        Input=Units(Input)
        Globals=GlobalConsts(kB=8.617333e-5,elecmass=1.0,eleccharge=1.0,ϵ0=8.065279e-9
        ,hbar=6.582e-1)
    elseif Input["Type"]=="2TM"
        Globals=GlobalConsts(kB=1.380649e-23,elecmass=9.10938e-31,eleccharge=1.60217663e-19
        ,ϵ0=8.854187e-12,hbar=1.054571e-34)
    end

    if Input["Dims"]==0
        Input["Length"]=1
    end

    simsettings=SimParam(ChemPot=Input["EF"],ElecPhon=Input["Gep"],
    ElecHeat=Input["Cel"],Dim=Input["Dims"],Output=Input["Output"],vers=Input["Type"],
    SimEnd=Input["SimEnd"])

    laser=LaserParam(FWHM=Input["FWHM"],ϕ=Input["Fluence"],delay=Input["LaserOff"],
    hv=Input["PhotEn"])

    matpat=MaterialParam(α=1/Input["ExtCof"],n=Input["AtomDens"],
    θ=Input["DebyeT"],ne=Input["Freeelec"],effmass=Input["Effmass"],FE=Input["FE"],
    ω=Input["Plasma"],DOSFile=Input["DOS"],γ=Input["gamma"],κ=Input["RTKappa"]
    ,L=Input["Length"],dz=Input["dz"],λ=0.0)

    
    slab=collect(0:matpat.dz:matpat.L)
    
    Simvariables=SimMaterialParam(Tel=fill(Input["Tel"],length(slab))
    ,Tph=fill(Input["Tel"],length(slab)),μ=Input["FE"],
    DOS=undef,lb=0.0,ub=0.0,NumElec=0.0,τ=0.0,g=Input["g"],zgrid=slab)

    return laser,matpat,simsettings,Globals,Simvariables

end

end
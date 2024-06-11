using ModelingToolkit
using ModelingToolkit: D_nounits as D, t_nounits as t

using DifferentialEquations,Plots,IfElse,Dierckx,DelimitedFiles,ForwardDiff # For testing

using .SymbolicsInterpolation

include("SimulationSetup.jl")


#Define types and strcts, each laser subtype is a functor for the ODE to build for the laser so new lasers require new subtypes
abstract type Laser end
abstract type LaserType end
@kwdef struct Gaussian <: LaserType 
    FWHM::Real
    Offset::Real
    Power::Real
    hv::Real
    Reflectivity::Real
    Transport::String
end
@kwdef struct Rectangular <: LaserType  
    FWHM::Real
    Offset::Real
    Power::Real
    hv::Real
    Reflectivity::Real
    Transport::String
end
@kwdef struct Secant <: LaserType
    FWHM::Real
    Offset::Real
    Power::Real
    hv::Real
    Reflectivity::Real
    Transport::String
end
@kwdef struct Lorentzian <: LaserType
    FWHM::Real
    Offset::Real
    Power::Real
    hv::Real
    Reflectivity::Real
    Transport::String
end

function define_laser_system(;lasertype::String,fwhm::Real,fluence::Real,
    photon_en::Real,offset::Real,reflectivity=0.0::Real,transport="Optical"::String)
    if lasertype == "Gaussian"
        return Gaussian(FWHM=fwhm,Offset=offset,Power=fluence,hv=photon_en,
        Reflectivity=reflectivity,Transport=transport)
    elseif lasertype == "Lorentzian"
        return Lorentzian(FWHM=fwhm,Offset=offset,Power=fluence,hv=photon_en,
        Reflectivity=reflectivity,Transport=transport)
    elseif lasertype == "Secant"
        return Secant(FWHM=fwhm,Offset=offset,Power=fluence,hv=photon_en,
        Reflectivity=reflectivity,Transport=transport)
    elseif lasertype == "Rectangular"
        return Rectangular(FWHM=fwhm,Offset=offset,Power=fluence,hv=photon_en,
        Reflectivity=reflectivity,Transport=transport)
    else
        print("The laser type you chose is currently implemented, the current options
        are Gaussian, Lorentzian, Secant and Rectangular")
    end
end

function define_laser_system(dict;lasertype=dict.laser::Laser_Type,fwhm=dict.FWHM::Real,
    fluence=dict.Fluence::Real,photon_en=Dict.hv::Real,offset=dict.Delay::Real,reflectivity=dict.Reflectivity::Real)

    if lasertype == "Gaussian"
        return Gaussian(FWHM=fwhm,Offset=offset,Power=fluence,hv=photon_en,
        Reflectivity=reflectivity)
    elseif lasertype == "Lorentzian"
        return Lorentzian(FWHM=fwhm,Offset=offset,Power=fluence,hv=photon_en,
        Reflectivity=reflectivity)
    elseif lasertype == "Secant"
        return Secant(FWHM=fwhm,Offset=offset,Power=fluence,hv=photon_en,
        Reflectivity=reflectivity)
    elseif lasertype == "Rectangular"
        return Rectangular(FWHM=fwhm,Offset=offset,Power=fluence,hv=photon_en,
        Reflectivity=reflectivity)
    else
        print("The laser type you chose is currently implemented, the current options
        are Gaussian, Lorentzian, Secant and Rectangular")
    end
end

function laser_factory(laser::LaserType,dims=Homogenous()::Dimension)
    temporal = laser()
    power = intensity()
    spatial = spatial_laser(laser,dims)
    return temporal*spatial*power
end

function intensity()
    @parameters R ϕ
    return (1-R)*ϕ
end

#The various temporal profiles
function (::Gaussian)()
    @parameters FWHM Offset    
    return sqrt(4*log(2)/pi)/FWHM*exp(-4*log(2)*(t-(2*FWHM)-Offset)^2/FWHM^2)
end

function (::Secant)()
    @parameters FWHM Offset   
    sec = sech(2*log(1+sqrt(2))*((t-(2*FWHM)-Offset)/FWHM))^2
   return log(1+sqrt(2))/FWHM * sec
end

function (::Lorentzian)()
    @parameters FWHM Offset
    lorent = (1+(4/(1+sqrt(2))*((t-(2*FWHM)-Offset)/FWHM)^2))^-2
    return 4*sqrt(sqrt(2)-1)/(pi*FWHM)*lorent
end

function (::Rectangular)()
    @parameters FWHM Offset
    return IfElse.ifelse(t≤Offset+4*FWHM,IfElse.ifelse(Offset≤t,1/(4*FWHM),0.0),0.0)
end

function spatial_laser(laser::LaserType,slab::Dimension)
    z_laser = spatial_z_laser(laser,slab)
    xy_laser = spatial_xy_laser(slab)
    return z_laser.*xy_laser
end

function spatial_z_laser(laser::LaserType,slab::Dimension)
    @parameters zgrid Z 
    if laser.Transport == "Ballistic"
        @parameters δb
        if typeof(slab) == Homogenous
            return 1/δb
        else
            return 1/(δb*(1-exp(-Z/δb))).*exp.(-zgrid./δb)
        end
    elseif laser.Transport == "Optical"
        @parameters ϵ
        if typeof(slab) == Homogenous
            return 1/ϵ
        else
            return 1/(ϵ*(1-exp(-Z/ϵ))).*exp.(-zgrid./ϵ)
        end
    elseif laser.Transport == "Combined"
        @parameters δb ϵ
        if typeof(slab) == Homogenous
            return 1/(δb+ϵ)
        else
            return 1/((δb+ϵ)*(1-exp(-Z/(δb+ϵ)))).*exp.(-zgrid./(δb+ϵ))
        end
    end
end

function spatial_xy_laser(slab::Dimension)
    @parameters X
    if typeof(slab) == Cylindrical
        @parameters xygrid
        return 1/(pi*X^2).*exp.(-xygrid.^2 ./X^2)
    elseif typeof(slab) == Cubic
        @parameters xgrid ygrid Y
        return 1/(pi^2*X^2*Y^2).*exp.((-xgrid.^2 ./X^2)+(-ygrid.^2 ./Y^2))
    else
        return 1
    end
end

function get_zgrid(slab)
    if slab.Dimension == 1
        return slab.spatial_grid
    elseif slab.Dimension == 2
        zvalues = [x[2] for x in slab.spatial_grid]
        return zvalues[:,1]
    elseif slab.Dimension == 3
        zvalues = [x[3] for x in slab.spatial_grid]
        return zvalues[:,1]
    end
end

function get_xygrid(slab)
    if slab.Dimension == 2
        xyvalues = [x[1] for x in slab.spatial_grid]
        return xyvalues[1,:]
    elseif slab.Dimension == 3
        xyvalues = [(x[1],x[2]) for x in slab.spatial_grid]
        return xyvalues[1,:]
    end
end
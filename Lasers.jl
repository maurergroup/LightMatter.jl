using ModelingToolkit, DynamicQuantities
using ModelingToolkit: D_nounits as D, t_nounits as t

using DifferentialEquations,Plots # For testing

#Define types and strcts, each laser subtype is a functor for the ODE to build for the laser so new lasers require new subtypes
abstract type Laser end
abstract type LaserType end
@kwdef struct Gaussian <: LaserType 
    FWHM::Real
    Offset::Real
    Power::Real
    hv::Real
    Reflectivity=0.0::Real
end
@kwdef struct Rectangular <: LaserType  
    FWHM::Real
    Offset::Real
    Power::Real
    hv::Real
    Reflectivity=0.0::Real
end
@kwdef struct Secant <: LaserType
    FWHM::Real
    Offset::Real
    Power::Real
    hv::Real
    Reflectivity=0.0::Real
end
@kwdef struct Lorentzian <: LaserType
    FWHM::Real
    Offset::Real
    Power::Real
    hv::Real
    Reflectivity=0.0::Real
end

struct laser_sim <: Laser
    laser_eq :: Num
    laser_params <:LaserType
end


#This is the function that is called by the user and handles all of the API to build the laser ODE
function LaserBuilder(lasertype<:Laser_Type;sim::simulation,mp::simulation,dims=Homogenous::Dimension)
    dim = dims(mp)
    laser_parameter=lasertype()
    
    return laser_sim(laser_eq,laser_params)
end

#The overearching laser equation which multiples the power by the spatial and temporal component
function Laser(laser_parameter)
    laser_parameter=
    P=Power(laser_parameter)
    z = Spatial(laser_parameter,)
    return P*z*tp
end

function Power(laser)
    return (1-laser.Reflectivity)*laser.Fluence
end

#The various temporal profiles
function (::Gaussian)(laser<:LaserType)    
    return sqrt(4*log(2)/pi)/laser.FWHM*exp(-4*log(2)*(t-(2*laser.FWHM)-laser.Offset)^2/laser.FWHM^2)
end

function (::Secant)(laser<:LaserType)
    sec = sech(2*log(1+sqrt(2))*((t-(2*laser.FWHM)-laser.Offset)/laser.FWHM))^2
   return log(1+sqrt(2))/laser.FWHM * sec
end

function (::Lorentzian)(laser<:LaserType)
    lorent = (1+(4/(1+sqrt(2))*((t-(2*laser.FWHM)-laser.Offset)/laser.FWHM)^2))^-2
    return 4*sqrt(sqrt(2)-1)/(pi*laser.FWHM)*lorent
end

function (::Rectangular)(laser<:LaserType)
    if laser.Offset ≤ t ≤ laser.Offset+(4*laser.FWHM)
        return 1/(4*laser.FWHM)
    else
        return 0.0
    end
end

function (::Homogenous)(matpat::MaterialParameters)
    return 1
end

function (::Linear)(matpat::MaterialParameters)
    return 
end
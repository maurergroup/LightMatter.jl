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
    photon_en::Real,offset::Real,reflectivity=0.0::Real)
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

function LaserBuilder(;laser::LaserType,sim::simulation,mp::simulation,dims=Homogenous::Dimension)
    temporal = laser(laser)
    power = Power(laser)
    if typeof(dims) == Homogenous
        return temporal*power
    else
        spatial = spatial_Laser(laser,dims,matpat)
        return temporal*spatial*power
    end
end

function Power(laser)
    return (1-laser.Reflectivity)*laser.Fluence
end

#The various temporal profiles
function (::Gaussian)(laser::LaserType)    
    return sqrt(4*log(2)/pi)/laser.FWHM*exp(-4*log(2)*(t-(2*laser.FWHM)-laser.Offset)^2/laser.FWHM^2)
end

function (::Secant)(laser::LaserType)
    sec = sech(2*log(1+sqrt(2))*((t-(2*laser.FWHM)-laser.Offset)/laser.FWHM))^2
   return log(1+sqrt(2))/laser.FWHM * sec
end

function (::Lorentzian)(laser::LaserType)
    lorent = (1+(4/(1+sqrt(2))*((t-(2*laser.FWHM)-laser.Offset)/laser.FWHM)^2))^-2
    return 4*sqrt(sqrt(2)-1)/(pi*laser.FWHM)*lorent
end

function (::Rectangular)(laser::LaserType)
    return IfElse.ifelse(laser.Offset≤t≤laser.Offset+4*laser.FWHM,1/(4*laser.FWHM),0.0)
end

function spatial_Laser(laser::LaserType,slab::Dimension,matpat::MaterialParameters)
    z_laser = spatial_z_Laser(laser::LaserType,slab::Dimension,matpat::MaterialParameters)
    xy_laser = spatial_xy_Laser(laser::LaserType,slab::Dimension,matpat::MaterialParameters)
    return z_laser.*xy_laser
end

function spatial_z_Laser(laser::LaserType,slab::Dimension,matpat::MaterialParameters)
    grid = get_zgrid(slab)
    if laser.Transport == "Ballistic"
        return 1/(matpat.ballistic*exp(-slab.Length/matpat.ballistic)).*exp.(-grid./matpat.ballistic)
    elseif laser.Transport == "Optical"
        return 1/(matpat.ballistic*exp(-slab.Length/matpat.ballistic)).*exp.(-grid./matpat.ballistic)
    elseif laser.Transport == "Combined"
        return 1/((matpat.ballistic+matpat.ϵ)*exp(-slab.Length/(matpat.ballistic+matpat.ϵ))).*exp.(-grid./(matpat.ballistic+matpat.ϵ))
    end
end

function spatial_xy_Laser(laser::LaserType,slab::Dimension,matpat::MaterialParameters)
    grid = get_xygrid(slab)
    if slab.Dimension == 2
        return 1/(pi*length(grid)^2).*exp.(-grid.^2 ./length(grid^2))
    elseif slab.Dimension == 3
        xgrid = [x[1] for x in grid]
        ygrid = [x[2] for x in grid]
        return 1/(pi*length(grid)^2).*exp.(-(xgrid.+ygrid).^2 ./length(grid^2))
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
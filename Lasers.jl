"""
    The supertype for the different laser types. Each laser type requires a struct which holds the laser parameters.
    Each laser type also has a functor defined for it which returns the equation describing how the laser evolves 
    over time. Lasers are centred on t = 0
"""
@kwdef struct Gaussian <: Laser
    FWHM::Real
    ϕ::Real
    hv::Real
    R::Real
    Transport::String
end
@kwdef struct Rectangular <: Laser  
    FWHM::Real
    ϕ::Real
    hv::Real
    R::Real
    Transport::String
end
@kwdef struct HyperbolicSecant <: Laser
    FWHM::Real
    ϕ::Real
    hv::Real
    R::Real
    Transport::String
end
@kwdef struct Lorentzian <: Laser
    FWHM::Real
    ϕ::Real
    hv::Real
    R::Real
    Transport::String
end
"""
    Outer constructor that determines which laser system to build. The parameters can be read in from InputFileControl.jl
    or you can manually declare them. Temporary: File support not supported
"""
function define_laser_system(Laser::Symbol;fwhm::Real,fluence::Real,
    photon_en::Real,reflectivity=0.0::Real,transport="Optical"::String)
    if Laser == :Gaussian
        return Gaussian(FWHM=fwhm,ϕ=fluence,hv=photon_en,
        R=reflectivity,Transport=transport)
    elseif Laser == :Lorentzian
        return Lorentzian(FWHM=fwhm,ϕ=fluence,hv=photon_en,
        R=reflectivity,Transport=transport)
    elseif Laser == :Secant
        return Secant(FWHM=fwhm,ϕ=fluence,hv=photon_en,
        R=reflectivity,Transport=transport)
    elseif Laser == :Rectangular
        return Rectangular(FWHM=fwhm,ϕ=fluence,hv=photon_en,
        R=reflectivity,Transport=transport)
    else
        print("The laser type you chose is currently not implemented, the current options
        are Gaussian, Lorentzian, Secant and Rectangular")
    end
end

function define_laser_system(dict;Laser=dict.laser::Symbol,fwhm=dict.FWHM::Real,
    fluence=dict.Fluence::Real,photon_en=Dict.hv::Real,reflectivity=dict.Reflectivity::Real)

    if Laser == "Gaussian"
        return Gaussian(FWHM=fwhm,ϕ=fluence,hv=photon_en,
        R=reflectivity)
    elseif Laser == "Lorentzian"
        return Lorentzian(FWHM=fwhm,ϕ=fluence,hv=photon_en,
        R=reflectivity)
    elseif Laser == "Secant"
        return HyperbolicSecant(FWHM=fwhm,ϕ=fluence,hv=photon_en,
        R=reflectivity)
    elseif Laser == "Rectangular"
        return Rectangular(FWHM=fwhm,ϕ=fluence,hv=photon_en,
        R=reflectivity)
    else
        print("The laser type you chose is currently implemented, the current options
        are Gaussian, Lorentzian, Secant and Rectangular")
    end
end
"""
    Factory builds the laser from three components defining the evolution in time(temporal), the power of the
    laser(power) and how the laser penetrates into a material(spatial). Returns a Num equation that holds the 
    relevant parameters and structure of the laser.
"""
function laser_factory(laser::Laser,mp::MaterialParameters,dims=Homogenous()::Dimension)
    temporal = laser(laser)
    power = intensity(laser)
    spatial = spatial_laser(laser,dims,mp)
    return temporal*spatial*power
end
"""
    Creates the parameters that represent the Reflectivity(R) of the material and the Fluence(ϕ) of the laser.
    Returns the component of the laser equation responsible for the laser power absorbed by the material.
"""
function intensity(laser)
    return (1-laser.R)*laser.ϕ
end
"""
    Normalised Gaussian temporal profile which creates the parameters for the Full-Width at Half-Maximum(FWHM).
"""
function (::Gaussian)(laser)
    return sqrt(4*log(2)/pi)/laser.FWHM*exp(-4*log(2)*t^2/laser.FWHM^2)
end
"""
    Normalised Hyperbolic Secant temporal profile which creates the parameters for the Full-Width at Half-Maximum(FWHM).
"""
function (::HyperbolicSecant)(laser)
    sec = sech(2*log(1+sqrt(2))*(t/laser.FWHM))^2
   return log(1+sqrt(2))/laser.FWHM * sec
end
"""
    Normalised Lorentzian temporal profile which creates the parameters for the Full-Width at Half-Maximum(FWHM).
"""
function (::Lorentzian)(laser)
    lorent = (1+(4/(1+sqrt(2))*(t/laser.FWHM)^2))^-2
    return 4*sqrt(sqrt(2)-1)/(pi*laser.FWHM)*lorent
end
"""
    Normalised Rectangular(continuous illumination) temporal profile which creates the parameters
    for the Full-Width at Half-Maximum(FWHM).
"""
function (::Rectangular)(laser)
    return IfElse.ifelse(t≤2*laser.FWHM,IfElse.ifelse(-2*laser.FWHM≤t,1/(4*laser.FWHM),0.0),0.0)
end
"""
    This builds the equation which defines how the laser energy is absorbed the further from
    the centre of the spot. This is split into a z and an x/y component.
"""
function spatial_laser(laser::Laser,slab::Dimension,mp::MaterialParameters)
    z_laser = spatial_z_laser(laser,slab,mp)
    xy_laser = spatial_xy_laser(slab)
    return z_laser.*xy_laser
end
"""
    The z-component of the laser profile is determined by the type of the slab object
    as well as whether ballistic or optical transport is being considered. It then creates
    the parameter for the ballistic length of electrons(δb) and the extinction coefficient (ϵ)
    depending on how the transport is modelled. It is assumed that the laser field is propogating
    along the z-axis
"""
function spatial_z_laser(laser::Laser,slab::Dimension,mp::MaterialParameters)
    @parameters zgrid Z 
    if laser.Transport == "Ballistic"
        if typeof(slab) == Homogenous
            return 1/mp.δb
        else
            return 1/(mp.δb*(1-exp(-Z/mp.δb))).*exp.(-zgrid./mp.δb)
        end
    elseif laser.Transport == "Optical"
        if typeof(slab) == Homogenous
            return 1/mp.ϵ
        else
            return 1/(mp.ϵ*(1-exp(-Z/mp.ϵ))).*exp.(-zgrid./mp.ϵ)
        end
    elseif laser.Transport == "Combined"
        if typeof(slab) == Homogenous
            return 1/(mp.ϵ+mp.ϵ)
        else
            return 1/((mp.ϵ+mp.ϵ)*(1-exp(-Z/(mp.ϵ+mp.ϵ)))).*exp.(-zgrid./(mp.ϵ+mp.ϵ))
        end
    end
end
"""
    Defines the equation for how the laser pulse diffuses perpendicular to it's incidence angle
    of incidence. Returns 1 unless the model is 2 or 3D in which it defines a parameter for the length
    of the slab defined by a capital letter and a parameter for the vector of slab positions in each
    dimension.
"""
function spatial_xy_laser(slab::Dimension)
    if typeof(slab) == Cylindrical
        @parameters xygrid R
        return 1/(pi*R^2).*exp.(-xygrid.^2 ./X^2)
    elseif typeof(slab) == Cubic
        @parameters xgrid ygrid Y X
        return 1/(pi^2*X^2*Y^2).*exp.((-xgrid.^2 ./X^2)+(-ygrid.^2 ./Y^2))
    else
        return 1
    end
end
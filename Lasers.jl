"""
    The supertype for the different laser types. Each laser type requires a struct which holds the laser parameters.
    Each laser type also has a functor defined for it which returns the equation describing how the laser evolves 
    over time. Lasers are centred on t = 0
"""
@kwdef struct Gaussian <: Laser
    FWHM::Real
    Power::Real
    hv::Real
    R::Real
    Transport::String
end
@kwdef struct Rectangular <: Laser  
    FWHM::Real
    Power::Real
    hv::Real
    R::Real
    Transport::String
end
@kwdef struct HyperbolicSecant <: Laser
    FWHM::Real
    Power::Real
    hv::Real
    R::Real
    Transport::String
end
@kwdef struct Lorentzian <: Laser
    FWHM::Real
    Power::Real
    hv::Real
    R::Real
    Transport::String
end
"""
    Outer constructor that determines which laser system to build. The parameters can be read in from InputFileControl.jl
    or you can manually declare them. Temporary: File support not supported
"""
function define_laser_system(Laser::Symbol;fwhm::Real,fluence::Real,
    photon_en::Real,R=0.0::Real,transport="Optical"::String)
    if Laser == :Gaussian
        return Gaussian(FWHM=fwhm,Power=fluence,hv=photon_en,
        R=R,Transport=transport)
    elseif Laser == :Lorentzian
        return Lorentzian(FWHM=fwhm,Power=fluence,hv=photon_en,
        R=R,Transport=transport)
    elseif Laser == :HyperbolicSecant
        return HyperbolicSecant(FWHM=fwhm,Power=fluence,hv=photon_en,
        R=R,Transport=transport)
    elseif Laser == :Rectangular
        return Rectangular(FWHM=fwhm,Power=fluence,hv=photon_en,
        R=R,Transport=transport)
    else
        print("The laser type you chose is currently not implemented, the current options
        are Gaussian, Lorentzian, HyperbolicSecant and Rectangular")
    end
end

#= function define_laser_system(dict;Laser=dict.laser::Symbol,fwhm=dict.FWHM::Real,
    fluence=dict.Fluence::Real,photon_en=Dict.hv::Real,R=dict.R::Real)
    if Laser == "Gaussian"
        return Gaussian(FWHM=fwhm,Power=fluence,hv=photon_en,
        R=R)
    elseif Laser == "Lorentzian"
        return Lorentzian(FWHM=fwhm,Power=fluence,hv=photon_en,
        R=R)
    elseif Laser == "HyperbolicSecant"
        return HyperbolicSecant(FWHM=fwhm,Power=fluence,hv=photon_en,
        R=R)
    elseif Laser == "Rectangular"
        return Rectangular(FWHM=fwhm,Power=fluence,hv=photon_en,
        R=R)
    else
        print("The laser type you chose is currently implemented, the current options
        are Gaussian, Lorentzian, HyperbolicSecant and Rectangular")
    end
end =#
"""
    Factory builds the laser from three components defining the evolution in time(temporal), the power of the
    laser(power) and how the laser penetrates into a material(spatial). Returns a Num equation that holds the 
    relevant parameters and structure of the laser.
"""
function laser_factory(las::Laser,dims=Homogeneous()::Dimension)
    temporal = las()
    power = intensity()
    spatial = spatial_laser(las,dims)
    return Expr(:call,:*,temporal,spatial,power)
end
"""
    Creates the parameters that represent the R(R) of the material and the Fluence(ϕ) of the laser.
    Returns the component of the laser equation responsible for the laser power absorbed by the material.
"""
function intensity()
    return :((1-las.R)*las.Power)
end
"""
    Normalised Gaussian temporal profile which creates the parameters for the Full-Width at Half-Maximum(FWHM).
"""
function (::Gaussian)()
    return :(sqrt(4*log(2)/pi)/las.FWHM*exp(-4*log(2)*t^2/las.FWHM^2))
end
"""
    Normalised Hyperbolic HyperbolicSecant temporal profile which creates the parameters for the Full-Width at Half-Maximum(FWHM).
"""
function (::HyperbolicSecant)()
   return :(log(1+sqrt(2))/las.FWHM * sech(2*log(1+sqrt(2))*(t/las.FWHM))^2)
end
"""
    Normalised Lorentzian temporal profile which creates the parameters for the Full-Width at Half-Maximum(FWHM).
"""
function (::Lorentzian)()
    lorent = :((1+(4/(1+sqrt(2))*(t/las.FWHM)^2))^-2)
    return :(4*sqrt(sqrt(2)-1)/(pi*las.FWHM)*$lorent)
end
"""
    Normalised Rectangular(continuous illumination) temporal profile which creates the parameters
    for the Full-Width at Half-Maximum(FWHM).
"""
function (::Rectangular)()
    return :(-2*las.FWHM≤t≤2*las.FWHM ? 1/(4*las.FWHM) : 0.0)
end

#
# To-Do Spatial Lasers
#
"""
    This builds the equation which defines how the laser energy is absorbed the further from
    the centre of the spot. This is split into a z and an x/y component.
"""
function spatial_laser(las::Laser,slab::Dimension)
    z_laser = spatial_z_laser(las,slab)
    return z_laser
end
"""
    The z-component of the laser profile is determined by the type of the slab object
    as well as whether ballistic or optical transport is being considered. It then creates
    the parameter for the ballistic length of electrons(δb) and the extinction coefficient (ϵ)
    depending on how the transport is modelled. It is assumed that the laser field is propogating
    along the z-axis
"""
function spatial_z_laser(las::Laser,slab::Dimension)
    las_vec = Vector{Expr}(undef,length(slab.grid))
    if las.Transport == "Ballistic"
        if typeof(slab) == Homogeneous
            return :(1/mp.δb)
        else
            zs = [x[end] for x in slab.grid]
            return :(1/(mp.δb*(1-exp(-length(slab.grid)/mp.δb))).*exp.(-$zs./mp.δb))
        end
    elseif las.Transport == "Optical"
        if typeof(slab) == Homogeneous
            return :(1/mp.ϵ)
        else
            l=slab.grid[end]
            z = :(dim.grid[i])
            las_vec= :(1/(mp.ϵ*(1-exp(-$l/mp.ϵ))).*exp.(-$z./mp.ϵ))
            return las_vec
        end
    elseif las.Transport == "Combined"
        if typeof(slab) == Homogeneous
            return :(1/(mp.δb+mp.ϵ))
        else
            zs = [x[end] for x in slab.grid]
            return :(1/((mp.δb+mp.ϵ)*(1-exp(-length(slab.grid)/(mp.δb+mp.ϵ)))).*exp.(-$zs[i]./(mp.δb+mp.ϵ)))
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
        return 1/(pi*R^2).*exp.(-xygrid.^2 ./X^2)
    elseif typeof(slab) == Cubic
        return 1/(pi^2*X^2*Y^2).*exp.((-xgrid.^2 ./X^2)+(-ygrid.^2 ./Y^2))
    else
        return 1
    end
end
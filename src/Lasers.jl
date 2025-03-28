@kwdef struct Gaussian <: Laser
    FWHM::Real
    Power::Real
    hv::Real
    Transport::String
end
@kwdef struct Rectangular <: Laser  
    FWHM::Real
    Power::Real
    hv::Real
    Transport::String
end
@kwdef struct HyperbolicSecant <: Laser
    FWHM::Real
    Power::Real
    hv::Real
    Transport::String
end
@kwdef struct Lorentzian <: Laser
    FWHM::Real
    Power::Real
    hv::Real
    Transport::String
end
"""
    define_laser_system(Laser::Symbol;fwhm::Real,fluence::Real,photon_en::Real,R=0.0::Real,transport="Optical"::String)
    The construcor which assembles the correct laser type for your problem. The currently implemented lasers are Gaussian, 
    Lorentzian, HyperbolicSecant and Rectangular. All lasers have the following structure,
    
    struct laser_type <:Laser
        FWHM::Real (The full-width half-maximum of the current laser, for rectangular it is half the duration either side of 0)
        Power::Real (The fluence of the laser before reflection)
        hv::Real (The frequency or photon energy of the laser)
        Transport::String (Either Optical, Ballistic or Combined for whether we are considering extinction coefficient, ballistic
                           length of electrons or both as ways that the laser energy reduces deeper into the slab.)
"""
function define_laser_system(Laser::Symbol;FWHM::Real,Power::Real,
    hv::Real,Transport="Optical"::String)
    if Laser == :Gaussian
        return Gaussian(FWHM=FWHM,Power=Power,hv=hv,Transport=Transport)
    elseif Laser == :Lorentzian
        return Lorentzian(FWHM=FWHM,Power=Power,hv=hv,Transport=Transport)
    elseif Laser == :HyperbolicSecant
        return HyperbolicSecant(FWHM=FWHM,Power=Power,hv=hv,Transport=Transport)
    elseif Laser == :Rectangular
        return Rectangular(FWHM=FWHM,Power=Power,hv=hv,Transport=Transport)
    else
        print("The laser type you chose is currently not implemented, the current options
        are Gaussian, Lorentzian, HyperbolicSecant and Rectangular")
    end
end
"""
    laser_factory(las::Laser,dims=Homogeneous()::Dimension)
    Factory builds the laser from three components defining the evolution in time(temporal), the power of the
    laser(power) and how the laser penetrates into a material(spatial). Returns an Expr that holds the 
    relevant parameters and structure of the laser.
"""
function laser_factory(las::Laser,dims=Homogeneous()::Dimension)
    temporal = las()
    power = :((1-mp.R)*las.Power)
    spatial = spatial_laser(las,dims)
    return Expr(:call,:*,temporal,spatial,power)
end
"""
    (::Gaussian)()
    Normalised Gaussian temporal profile
"""
function (::Gaussian)()
    return :(sqrt(4*log(2)/pi)/las.FWHM*exp(-4*log(2)*t^2/las.FWHM^2))
end
"""
    (::HyperbolicSecant)()
    Normalised Hyperbolic Secant temporal profile
"""
function (::HyperbolicSecant)()
   return :(log(1+sqrt(2))/las.FWHM * sech(2*log(1+sqrt(2))*(t/las.FWHM))^2)
end
"""
    (::Lorentzian)()
    Normalised Lorentzian temporal profile
"""
function (::Lorentzian)()
    lorent = :((1+(4/(1+sqrt(2))*(t/las.FWHM)^2))^-2)
    return :(4*sqrt(sqrt(2)-1)/(pi*las.FWHM)*$lorent)
end
"""
    (::Rectangular)()
    Normalised Rectangular temporal profile with illumination lasting 2*FWHM either side of 0 fs.
"""
function (::Rectangular)()
    return :(-2*las.FWHM≤t≤2*las.FWHM ? 1/(4*las.FWHM) : 0.0)
end
"""
    spatial_laser(las::Laser,slab::Dimension)
    This builds the equation which defines how the laser energy is absorbed the further from
    the centre of the spot. This is currently only implemented in the z direction
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
            return :(1/(mp.ϵ*(1-exp(-$l/mp.ϵ)))*exp.(-dim.grid[i]./mp.ϵ))
        end
    elseif las.Transport == "Combined"
        if typeof(slab) == Homogeneous
            return :(1/(mp.δb+mp.ϵ))
        else
            l=slab.grid[end]
            return :(1/((mp.δb+mp.ϵ)*(1-exp(-$l/(mp.δb+mp.ϵ)))).*exp.(-dim.grid[z]./(mp.δb+mp.ϵ)))
        end
    end
end
"""
    spatial_xy_laser(slab::Dimension)
    To-Do : Not currently working and implemented yet. 
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
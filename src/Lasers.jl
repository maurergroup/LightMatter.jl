"""
    laser_factory(las::Laser,dims=Homogeneous()::Dimension)
    Factory builds the laser from three components defining the evolution in time(temporal), the power of the
    laser(power) and how the laser penetrates into a material(spatial). Returns an Expr that holds the 
    relevant parameters and structure of the laser.
"""
function laser_factory(sim::Simulation)
    temporal = temporal_laser(sim)
    power = :((1-sim.laser.R)*sim.laser.ϕ)
    spatial = spatial_laser(sim)
    return Expr(:call,:*,temporal,spatial,power)
end
"""
    temporal_laser(sim::Simulation)
    Returns the expression for the requested laser envelope function. Current options are
    :Gaussian, :HyperbolicSecant, :Lorentzian and :Rectangular.
"""
function temporal_laser(sim::Simulation)
    type = sim.laser.envelope
    if type == :Gaussian
        return :(sqrt(4*log(2)/pi)/las.FWHM*exp(-4*log(2)*t^2/las.FWHM^2))
    elseif type == :HyperbolicSecant
        return :(log(1+sqrt(2))/las.FWHM * sech(2*log(1+sqrt(2))*(t/las.FWHM))^2)
    elseif type == :Lorentzian
        lorent = :((1+(4/(1+sqrt(2))*(t/las.FWHM)^2))^-2)
        return :(4*sqrt(sqrt(2)-1)/(pi*las.FWHM)*$lorent)
    elseif type == :Rectangular
        return :(-2*las.FWHM≤t≤2*las.FWHM ? 1/(4*las.FWHM) : 0.0)
    end
end
"""
    spatial_laser(las::Laser,slab::Dimension)
    This builds the equation which defines how the laser energy is absorbed the further from
    the centre of the spot. This is currently only implemented in the z direction
"""
function spatial_laser(sim::Simulation)
    z_laser = spatial_z_laser(sim)
    return z_laser
end
"""
    The z-component of the laser profile is determined by the type of the slab object
    as well as whether ballistic or optical transport is being considered. It then creates
    the parameter for the ballistic length of electrons(δb) and the extinction coefficient (ϵ)
    depending on how the transport is modelled. It is assumed that the laser field is propogating
    along the z-axis
"""
function spatial_z_laser(sim::Simulation)
    if sim.laser.Transport == "Ballistic"
        if sim.structure.dimension.length == 1
            return :(1/sim.laser.δb)
        else
            l = sim.structure.dimension.grid[end]
            return :(1/(sim.laser.δb*(1-exp(-$l/sim.laser.δb)))*exp.(-sim.structure.dimension.grid[i]./sim.laser.δb))
        end
    elseif sim.laser.Transport == "Optical"
        if sim.structure.dimension.length == 1
            return :(1/sim.laser.ϵ)
        else
            l = sim.structure.dimension.grid[end]
            return :(1/(sim.laser.ϵ*(1-exp(-$l/sim.laser.ϵ)))*exp.(-sim.structure.dimension.grid[i]./sim.laser.ϵ))
        end
    elseif sim.laser.Transport == "Combined"
        if sim.structure.dimension.length == 1
            return :(1/(sim.laser.δb+sim.laser.ϵ))
        else
            l = sim.structure.dimension.grid[end]
            return :(1/((sim.laser.δb+sim.laser.ϵ)*(1-exp(-$l/(sim.laser.δb+sim.laser.ϵ)))).*exp.(-sim.structure.dimension.grid[i]./(sim.laser.δb+sim.laser.ϵ)))
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
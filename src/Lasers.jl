"""
    laser_factory(sim::Simulation)
    
    Assembles the user desired laser expression

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Expression for the power of the laser as a function of time and depth
"""
function laser_factory(sim::Simulation)
    temporal = temporal_laser(sim)
    power = :((1-sim.laser.R) * sim.laser.ϕ)
    spatial = spatial_laser(sim)
    return Expr(:call, :*, temporal, spatial, power)
end
"""
    temporal_laser(sim::Simulation)
    
    Assembles the expression for the user desired temporal shape of the laser
    Currently implemented:
    - :Gaussian : Gaussian laser shape
    - :Lorentzian : Lorentzian laser shape
    - :HyperbolicSecant : Hyperbolic secant laser shape
    - :Rectangular : Constant illumination style method. The laser is on for 2*FWHM at either side of 0.0 fs

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Expression for the temporal shape of the laser
"""
function temporal_laser(sim::Simulation)
    type = sim.laser.envelope
    if type == :Gaussian
        return :(sqrt(4*log(2)/pi)/sim.laser.FWHM*exp(-4*log(2)*t^2/sim.laser.FWHM^2))
    elseif type == :HyperbolicSecant
        return :(log(1+sqrt(2))/sim.laser.FWHM * sech(2*log(1+sqrt(2))*(t/sim.laser.FWHM))^2)
    elseif type == :Lorentzian
        lorent = :((1+(4/(1+sqrt(2))*(t/sim.laser.FWHM)^2))^-2)
        return :(4*sqrt(sqrt(2)-1)/(pi*sim.laser.FWHM)*$lorent)
    elseif type == :Rectangular
        return :(-2*sim.laser.FWHM≤t≤2*sim.laser.FWHM ? 1/(4*sim.laser.FWHM) : 0.0)
    end
end
"""
    spatial_laser(sim::Simulation)
    
    Assembles the expression for the laser penetration depth. Currently only works in the z-axis

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Expression for the spatial shape of the laser
"""
function spatial_laser(sim::Simulation)
    z_laser = spatial_z_laser(sim)
    return z_laser
end
"""
    spatial_z_laser(sim::Simulation)
    
    Assembles the expression for the laser penetration depth. Currently only works in the z-axis

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Expression for the spatial shape of the laser
"""
function spatial_z_laser(sim::Simulation)
    if sim.laser.Transport == :ballistic
        if sim.structure.dimension.length == 1
            return :(1/sim.laser.δb)
        else
            l = sim.structure.dimension.grid[end]
            return :(1/(sim.laser.δb*(1-exp(-$l/sim.laser.δb)))*exp.(-sim.structure.dimension.grid[i]./sim.laser.δb))
        end
    elseif sim.laser.Transport == :optical
        if sim.structure.dimension.length == 1
            return :(1/sim.laser.ϵ)
        else
            l = sim.structure.dimension.grid[end]
            return :(1/(sim.laser.ϵ*(1-exp(-$l/sim.laser.ϵ)))*exp.(-sim.structure.dimension.grid[i]./sim.laser.ϵ))
        end
    elseif sim.laser.Transport == :combined
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
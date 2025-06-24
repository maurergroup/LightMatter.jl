"""
    laser_factory(sim::Simulation)
    
    Assembles the user desired laser expression

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Expression for the power of the laser as a function of time and depth
"""
function laser_factory(sim::Simulation)
    temporal = temporal_laser(sim.laser)
    power = :((1-sim.laser.R) * sim.laser.ϕ)
    spatial = spatial_laser(sim)
    return Expr(:call, :*, temporal, spatial, power)
end
"""
    temporal_laser(las::Laser)
    
    Assembles the expression for the user desired temporal shape of the laser
    Currently implemented:
    - :Gaussian : Gaussian laser shape
    - :Lorentzian : Lorentzian laser shape
    - :HyperbolicSecant : Hyperbolic secant laser shape
    - :Rectangular : Constant illumination style method. The laser is on for 2*FWHM at either side of 0.0 fs

    # Arguments
    - 'las': Settings and parameters of the laser

    # Returns
    - Expression for the temporal shape of the laser
"""
function temporal_laser(laser::Laser)
    type = laser.envelope
    if type == :Gaussian
        return :(sqrt(4*log(2)/pi) / sim.laser.FWHM * exp(-4*log(2)*t^2/sim.laser.FWHM^2))
    elseif type == :HyperbolicSecant
        return :(log(1+sqrt(2)) / sim.laser.FWHM * sech(2*log(1+sqrt(2))*(t/sim.laser.FWHM))^2)
    elseif type == :Lorentzian
        lorent = :((1+(4/(1+sqrt(2)) * (t/sim.laser.FWHM)^2))^-2)
        return :(4 * sqrt(sqrt(2)-1) / (pi*sim.laser.FWHM) * $lorent)
    elseif type == :Rectangular
        return :(-2*sim.laser.FWHM ≤ t ≤ 2*sim.laser.FWHM ? 1/(4*sim.laser.FWHM) : 0.0)
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
    
    Defines how the laser energy changes with depth into the sample
    Currently implemented:
    - :ballistic : Ballistic depth of electrons defines spatial change, controlled by δb
    - :optical : Absorbtion depth of the laser controls the spatial change, controlled by ϵ (1/α)
    - :combined : Both ballistic and optical parameters enabled, sums the two effects together

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Expression for the spatial shape of the laser
"""
function spatial_z_laser(sim::Simulation)
    if sim.laser.Transport == :ballistic
        if sim.structure.dimension.length == 1
            return :(1 ./ sim.laser.δb)
        else
            l = sim.structure.dimension.grid[end]
            return :(1 ./ (sim.laser.δb.*(1 .-exp.(-$l./sim.laser.δb))) .* exp.(-sim.structure.dimension.grid[i]./sim.laser.δb))
        end
    elseif sim.laser.Transport == :optical
        if sim.structure.dimension.length == 1
            return :(1 ./ sim.laser.ϵ)
        else
            l = sim.structure.dimension.grid[end]
            return :(1 ./(sim.laser.ϵ .* (1 .-exp.(-$l./sim.laser.ϵ))) .* exp.(-sim.structure.dimension.grid[i]./sim.laser.ϵ))
        end
    elseif sim.laser.Transport == :combined
        if sim.structure.dimension.length == 1
            return :(1 ./ (sim.laser.δb .+ sim.laser.ϵ))
        else
            l = sim.structure.dimension.grid[end]
            return :(1 ./ ((sim.laser.δb.+sim.laser.ϵ) * (1 .-exp.(-$l./(sim.laser.δb.+sim.laser.ϵ)))) .* exp.(-sim.structure.dimension.grid[i]./(sim.laser.δb.+sim.laser.ϵ)))
        end
    end
end
"""
    WIP: spatial_xy_laser(sim::Simulation)
    
    Defines how the laser energy changes with distance from the spot centre. Do not use as simulations
    aren't constructed for more than 1D.

    Currently implemented:
    - :Cylindrical : Circular shaped sample
    - :Cubic : Cubic shaped sample

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Expression for the spatial shape of the laser in the xy plane
"""
function spatial_xy_laser(sim::Simulation)
    if typeof(slab) == Cylindrical
        return 1/(pi*R^2).*exp.(-xygrid.^2 ./X^2)
    elseif typeof(slab) == Cubic
        return 1/(pi^2*X^2*Y^2).*exp.((-xgrid.^2 ./X^2)+(-ygrid.^2 ./Y^2))
    else
        return 1
    end
end

function calculate_laser_fields(las::Laser)
    return get_laser_fields(las)
end

function get_laser_fields(las)
    if las.hv isa Matrix
        power = :(sim.laser.hv[:,2]*abs(sim.laser.hv[2,1]-sim.laser.hv[1,1]))  
    else
        power = :(sim.laser.hv)
    end
    if las.envelope == :Rectangular
        E_0 = :(-2*sim.laser.FWHM ≤ t ≤ 2*sim.laser.FWHM ? sqrt.(2*sim.laser.ϕ.*$power ./ (Constants.c*Constants.ϵ0*4*sim.laser.FWHM*sim.laser.n)) : 0.0)
    elseif las.envelope == :Gaussian
        E_0 = :(sqrt.(2*sim.laser.ϕ*sqrt(4*log(2)).*$power./(Constants.c*Constants.ϵ0*sim.laser.FWHM*sim.laser.n*sqrt(pi))).*exp(-2*log(2)*t^2/sim.laser.FWHM^2))
    end

    if las.hv isa Matrix
        E = :(sum($E_0 .* cos.(sim.laser.hv * 1/(Constants.ħ*2*pi) * t)))
    else
        E = :($E_0 .* cos.(sim.laser.hv * 1/(Constants.ħ*2*pi) * t))
    end
    B = :($E ./ Constants.c)
    return Fields(E, B)
end
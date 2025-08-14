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
    power = :((1 .-sim.laser.R) .* sim.laser.ϕ)
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
    if ndims(sim.structure.dimension.grid) == 3
        xylaser = spatial_xy_laser(sim) 
        return :(z_laser * xylaser)
    else
        return z_laser
    end
end

function heaviside(t)
   return 0.5 * (sign(t) + 1)
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
#= function spatial_z_laser(sim::Simulation)
    decay = spatial_laser_decay(sim)
    if sim.structure.dimension.length == 1
        return :(1 ./ $decay)
    else
        return :(1 ./ $decay .* exp.(-sim.structure.dimension.grid[i]./$decay))
    end
end =#

function spatial_z_laser(sim::Simulation)
    decay = spatial_laser_decay(sim)
    const_expr = :(1 ./ $decay)
    if sim.structure.Elemental_System != 1
        exprs = antennareactor_laserdecay(sim)
        return Expr(:call, :+, exprs...)
    else
        return :($const_expr .* exp.(-sim.structure.dimension.grid[i]./$decay))
    end 
end

function antennareactor_laserdecay(sim::Simulation)
    depths = vcat(sim.structure.dimension.InterfaceHeight, sim.structure.dimension.grid[end])

    layer_exprs = Expr[]
    for i in 1:sim.structure.Elemental_System
        ϵ = sim.laser.ϵ[i]
        z0 = (i == 1) ? 0.0 : depths[i-1]   # start of layer
        zend = depths[i]                     # end of layer

        # attenuation from all previous layers
        atten_prev = exp(-sum(diff([0.0; depths[1:i-1]] ./ sim.laser.ϵ[1:i-1])))
        
        expr = :( (1 / $ϵ) * $atten_prev * exp(-(sim.structure.dimension.grid[i] - $z0) / $ϵ) *
                  (heaviside(sim.structure.dimension.grid[i] - $z0) * heaviside($zend - sim.structure.dimension.grid[i])) )
        push!(layer_exprs, expr)
    end
    
    return layer_exprs
end

function spatial_laser_decay(sim::Simulation)
    if sim.laser.Transport == :ballistic
        return :(sim.laser.δb)
    elseif sim.laser.Transport == :optical
        return :(sim.laser.ϵ)
    elseif sim.laser.Transport == :combined
        return  :(sim.laser.δb .+ sim.laser.ϵ)
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
    if sim.structure.dimension.grid
        return 1/(pi*R^2).*exp.(-xygrid.^2 ./X^2)
    elseif typeof(slab) == Cubic
        return 1/(pi^2*X^2*Y^2).*exp.((-xgrid.^2 ./X^2)+(-ygrid.^2 ./Y^2))
    else
        return 1
    end
end
"""
    get_laser_fields(las::Laser)
    
    Returns an expression for the electric field induced by the laser
    Currently implemented:
    - :Gaussian
    - :Rectangular

    # Arguments
    - 'las': Laser settings 

    # Returns
    - Expression for the electric field of the laser
"""
function get_laser_fields(las)
    if las !== nothing
        if las.hv isa Matrix
            power = :(sim.laser.hv[:,2])  
        else
            power = :(1.0)
        end
        if las.envelope == :Rectangular
            E_0 = :(-2*sim.laser.FWHM ≤ t ≤ 2*sim.laser.FWHM ? sqrt.(2*sim.laser.ϕ.*$power ./ (Constants.c*Constants.ϵ0*4*sim.laser.FWHM*sim.laser.n)) : 0.0)
        elseif las.envelope == :Gaussian
            E_0 = :(sqrt.(2*sim.laser.ϕ*sqrt(4*log(2)).*$power./(Constants.c*Constants.ϵ0*sim.laser.FWHM*sim.laser.n*sqrt(pi))).*exp(-2*log(2)*t^2/sim.laser.FWHM^2))
        end

        if las.hv isa Matrix
            E = [:(sum($E_0 .* cos.(sim.laser.hv * 1/(Constants.ħ*2*pi) * t))), :(0.0+0.0), :(0.0+0.0)]
        else
            E = [:($E_0 .* cos.(sim.laser.hv * 1/(Constants.ħ*2*pi) * t)), :(0.0+0.0), :(0.0+0.0)]
        end
        B = [:(0.0+0.0), :($E ./ Constants.c), :(0.0+0.0)]
        return Fields(E, B)
    else 
        return Fields(:(0.0), :(0.0))
    end
end
"""
    photon_energytofrequency(hv::Real)
    
    # Arguments
    - 'hv': Photon energy in eV

    # Returns
    - Photon frequency in 1 / fs 
"""
function photon_energytofrequency(hv)
    return hv / Constants.ħ
end
"""
    E_magnitude(las_field::Vector{Expr}, ext_field::Vector{Expr})
    
    # Arguments
    - 'las_field': Electric field vector from the laser
    - 'ext_field': Electric field vector from external sources

    # Returns
    - Expression for the total magntitude of electric fields
"""
function E_magnitude(las_field, ext_field)
    sum = :(las_field .+ ext_field)
    return :(sqrt($(sum[1])^2 + $(sum[2])^2 + $(sum[3])^2))
end
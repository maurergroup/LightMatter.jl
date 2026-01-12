using Pkg
cd(raw"C:\Users\u5522838\Dropbox\PhD\Light-Matter Interactions")
Pkg.activate(".")
using Unitful, CairoMakie #Temporary


abstract type Geometry end

function build_Geometry(geometry, length, grid, spacing, interface)
    if geometry == :OneDGeometry
        return build_OneDGeometry(grid, interface)
    elseif geometry == :CoreShell
        return constructor_2D(CoreShell, length, spacing, interface)
    elseif geometry == :CoreSatellite
        return constructor_2D(CoreSatellite, length, spacing, interface)
    elseif geometry == :AntennaReactor
        return constructor_2D(AntennaReactor, length, spacing, interface)
    else
        error("Geometry type not recognized.")
    end
end

@kwdef struct OneDGeometry <: Geometry
    length::Int # The length of the grid, not the depth of the slab
    grid::AbstractArray{Float64} # The grid the simulation is solved over
    spacing::Union{Float64, Vector{Float64}} #The spacing between grid points
    InterfaceHeight::Union{Float64, Vector{Float64}} # Height sorted list of the interfaces between materials
end
"""
    build_OneDGeometry(grid=[0.0]::AbstractArray{Float64}, cutoff=0.0::Union{Vector{Float64},Float64})

    Outer constructor function to assemble the OneDGeometry struct. The user provides an evenly spaced grid 
    and sorted list of interface heights for antenna-reactor complexes. The user must ensure the length of 
    cutoff = Elemental_System - 1. No unit conversion is performed when assembling this struct.
    Defaults allow any unneccessary parameters for users simulation to be ignored.

    # Arguments
    - 'grid': unit = nm: Vector representing spatial grid. If [0.0] then homogenous (0D) calculation is performed
    - 'cutoff': unit = nm:  Sorted list of all interface heights. Only used when Elemental_System > 1. 

    # Returns
    - The OneDGeometry struct with the users grid and interface heights
"""
function build_OneDGeometry(depth, spacing, InterfaceHeight::Union{Vector{Float64},Float64} = 0.0)
    grid = collect(0.0 : convert_units(u"nm", spacing) : convert_units(u"nm", depth))
    L = length(grid)
    return OneDGeometry(L, grid, spacing, InterfaceHeight)
end

abstract type TwoDGeometry <: Geometry end

@kwdef struct CoreShell <: TwoDGeometry
    grid::Array{Tuple{Float64,Float64},2}  # (r, θ) in polar coordinates
    mask::Array{Int,2} # Determines the material at each grid point
    lengths::Tuple{Int, Int}  # (n_radial, n_angular)
    InterfaceRadius::Float64 # Distance from centre of sphere to edge of core material
end
@kwdef struct CoreSatellite <: TwoDGeometry
    grid::Array{Tuple{Float64,Float64},2}  # (x, y) in Cartesian coordinates
    mask::Array{Int,2} # Determines the material at each grid point
    lengths::Tuple{Int, Int}  # (n_x, n_y)
    core_radius::Float64
    satellite_radius::Float64
    n_satellites::Int
end

@kwdef struct AntennaReactor <: TwoDGeometry
    grid::Array{Tuple{Float64,Float64},2}  # (x, z) in Cartesian coordinates (z = height)
    mask::Array{Int,2} # Determines the material at each grid point
    lengths::Tuple{Int, Int}  # (n_x, n_z)
    reactor_radius::Float64
    antenna_height::Float64
end

function build_mask(grid, interface)
    return [searchsortedlast(interface, x) + 1 for x in grid]
end

function build_CoreShell(radius, spacing, InterfaceRadius::Union{Vector{Float64},Float64} = 0.0)
    grid = radial_meshgrid(radius::Real, spacing::Real)
    mask = coreshell_mask(grid, InterfaceRadius)
    return CoreShell(grid, mask, size(grid), InterfaceRadius)
end

function radial_meshgrid(radius::Real, spacing::Real)
    # Create radial grid (r from 0 to radius)
    r_coords = 0.0:spacing:radius*2
    # Create angular grid (θ from 0 to 2π)
    # Number of angular points scales with radius to maintain reasonable resolution
    n_theta = max(4, Int(ceil(2 * π * radius / spacing)))
    theta_coords = collect(0.0:(2*π)/n_theta:(2*π - (2*π)/n_theta))
    
    # Create 2D polar grid
    grid = [(r, θ) for r in r_coords, θ in theta_coords]
    return grid
end

function coreshell_mask(grid, interface::Real)
    # Get maximum radius from grid
    max_r = maximum(r for (r, θ) in grid)
    # Core radius is inner circle: max_r - interface (shell thickness)
    core_radius = max_r - interface
    mask = zeros(Int, size(grid))
    
    for j in axes(grid, 2), i in axes(grid, 1)
        r = grid[i, j][1]  # Radial coordinate
        
        if r <= core_radius
            mask[i, j] = 2   # Core region (inner circle)
        else
            mask[i, j] = 1   # Shell region (outer ring, thickness = interface)
        end
    end
    return mask
end

function plot_coreshell_circles(coreshell::CoreShell)
    """
    Plots the core-shell geometry as concentric circles in Cartesian coordinates
    """
    max_r = maximum(r for (r, θ) in coreshell.grid)
    core_radius = max_r - coreshell.InterfaceRadius
    
    fig = Figure(size = (600, 600))
    ax = Axis(fig[1, 1], aspect = DataAspect())
    
    # Create circle points
    θ_circle = collect(0:0.01:2π)
    
    # Plot core circle
    core_x = core_radius .* cos.(θ_circle)
    core_y = core_radius .* sin.(θ_circle)
    lines!(ax, core_x, core_y, color = :red, linewidth = 2, label = "Core boundary")
    
    # Plot outer shell circle
    shell_x = max_r .* cos.(θ_circle)
    shell_y = max_r .* sin.(θ_circle)
    lines!(ax, shell_x, shell_y, color = :blue, linewidth = 2, label = "Shell boundary")
    
    # Fill regions with transparency
    poly!(ax, Point2f.(zip(core_x, core_y)), color = (:red, 0.2), label = "Core")
    poly!(ax, Point2f.(zip(shell_x, shell_y)), color = (:blue, 0.2), label = "Shell")
    
    ax.xlabel = "x (nm)"
    ax.ylabel = "y (nm)"
    ax.title = "Core-Shell Structure"
    xlims!(ax, -max_r*1.1, max_r*1.1)
    ylims!(ax, -max_r*1.1, max_r*1.1)
    axislegend(ax, position = :rt)
    
    return fig
end

function build_CoreSatellite(core_radius::Real, satellite_radius::Real, spacing::Real)
    """
    Builds a CoreSatellite geometry with a central core and 4 satellites at Cartesian poles.
    Only creates half the structure for performance (x >= 0).
    """
    # Grid range: core + 2 satellite radii on each side
    # Only half structure: x from 0 to core + 2*satellite_radius, y from -(core + satellite_radius) to (core + satellite_radius)
    range_x = core_radius + 2*satellite_radius
    range_y = core_radius + satellite_radius
    x_coords = 0.0:spacing:range_x
    y_coords = -range_y:spacing:range_y
    
    # Create 2D Cartesian grid (half structure)
    grid = [(x, y) for x in x_coords, y in y_coords]
    mask = coresatellite_mask(grid, core_radius, satellite_radius)
    
    return CoreSatellite(grid, mask, size(grid), core_radius, satellite_radius, 4)
end


function coresatellite_mask(grid::Array{Tuple{Float64,Float64},2}, core_radius::Real, satellite_radius::Real)
    """
    Creates mask for core-satellite geometry with 4 satellites at Cartesian poles.
    Label 2 = core, Label 1 = satellites, Label 0 = void
    Satellites positioned at: (+core_radius, 0), (-core_radius, 0), (0, +core_radius), (0, -core_radius)
    """
    mask = zeros(Int, size(grid))
    
    # Define satellite positions at Cartesian poles
    satellite_positions = [
        (core_radius, 0.0),      # +x pole
        (-core_radius, 0.0),     # -x pole
        (0.0, core_radius),      # +y pole
        (0.0, -core_radius)      # -y pole
    ]
    
    for j in axes(grid, 2), i in axes(grid, 1)
        x, y = grid[i, j]
        r_center = sqrt(x^2 + y^2)
        
        # Check if in core
        if r_center <= core_radius
            mask[i, j] = 2
        else
            # Check if in any satellite
            for (sat_x, sat_y) in satellite_positions
                dist_to_sat = sqrt((x - sat_x)^2 + (y - sat_y)^2)
                
                if dist_to_sat <= satellite_radius
                    mask[i, j] = 1
                    break
                end
            end
        end
    end
    return mask
end

function plot_coresatellite_structure(cs::CoreSatellite)
    """
    Plots the core-satellite geometry with core and 4 satellites at Cartesian poles.
    Shows only the half structure (x >= 0) for performance.
    """
    fig = Figure(size = (600, 600))
    ax = Axis(fig[1, 1], aspect = DataAspect())
    
    # Create circle points
    θ_circle = collect(0:0.01:2π)
    
    # Plot core circle (full for reference, but only half will be visible on grid)
    core_x = cs.core_radius .* cos.(θ_circle)
    core_y = cs.core_radius .* sin.(θ_circle)
    lines!(ax, core_x, core_y, color = :red, linewidth = 2, label = "Core")
    poly!(ax, Point2f.(zip(core_x, core_y)), color = (:red, 0.3))
    
    # Plot 4 satellites at Cartesian poles
    satellite_positions = [
        (cs.core_radius, 0.0),      # +x pole
        (-cs.core_radius, 0.0),     # -x pole
        (0.0, cs.core_radius),      # +y pole
        (0.0, -cs.core_radius)      # -y pole
    ]
    
    for (sat_x, sat_y) in satellite_positions
        # Draw satellite circle
        sat_circle_x = sat_x .+ cs.satellite_radius .* cos.(θ_circle)
        sat_circle_y = sat_y .+ cs.satellite_radius .* sin.(θ_circle)
        lines!(ax, sat_circle_x, sat_circle_y, color = :blue, linewidth = 2)
        poly!(ax, Point2f.(zip(sat_circle_x, sat_circle_y)), color = (:blue, 0.3))
    end
    
    ax.xlabel = "x (nm)"
    ax.ylabel = "y (nm)"
    ax.title = "Core-Satellite Structure (Half, Poles)"
    
    # Set limits for half structure view (x >= 0)
    limit_x = cs.core_radius + 2*cs.satellite_radius
    limit_y = cs.core_radius + cs.satellite_radius
    xlims!(ax, -0.5*limit_x, limit_x*1.1)  # Show small portion of -x for context
    ylims!(ax, -limit_y*1.1, limit_y*1.1)
    axislegend(ax, position = :rt)
    
    return fig
end

function build_AntennaReactor(reactor_radius::Real, antenna_height::Real, spacing::Real)
    """
    Builds an AntennaReactor geometry with a reactor base and antenna on top (vertical configuration).
    """
    # Grid range: reactor diameter in x, reactor + antenna in z
    x_coords = -reactor_radius:spacing:reactor_radius
    z_coords = -reactor_radius:spacing:antenna_height
    
    # Create 2D Cartesian grid (x, z)
    grid = [(x, z) for x in x_coords, z in z_coords]
    mask = antennareactor_mask(grid, reactor_radius, antenna_height)
    
    return AntennaReactor(grid, mask, size(grid), reactor_radius, antenna_height)
end


function antennareactor_mask(grid::Array{Tuple{Float64,Float64},2}, reactor_radius::Real, antenna_height::Real)
    """
    Creates mask for antenna-reactor geometry.
    Label 2 = reactor base, Label 1 = antenna nanoparticle, Label 0 = void
    """
    mask = zeros(Int, size(grid))
    antenna_radius = reactor_radius * 0.3  # Antenna is 30% of reactor radius
    antenna_center_z = antenna_height / 2.0
    
    for j in axes(grid, 2), i in axes(grid, 1)
        x, z = grid[i, j]
        
        # Reactor: semicircle below z=0
        if z < 0 && sqrt(x^2 + z^2) <= reactor_radius
            mask[i, j] = 2
        # Antenna: spherical cap above z=0
        elseif z > 0 && z <= antenna_height
            dist_to_antenna = sqrt(x^2 + (z - antenna_center_z)^2)
            if dist_to_antenna <= antenna_radius
                mask[i, j] = 1
            end
        end
    end
    return mask
end

function plot_antennareactor_structure(ar::AntennaReactor)
    """
    Plots the antenna-reactor geometry with reactor base and antenna on top.
    """
    fig = Figure(size = (600, 600))
    ax = Axis(fig[1, 1], aspect = DataAspect())
    
    # Create semi-circle points for reactor (below z=0)
    θ_reactor = collect(π:0.01:2π)  # Lower semicircle
    reactor_x = ar.reactor_radius .* cos.(θ_reactor)
    reactor_z = ar.reactor_radius .* sin.(θ_reactor)
    lines!(ax, reactor_x, reactor_z, color = :red, linewidth = 2, label = "Reactor")
    poly!(ax, Point2f.(zip(reactor_x, reactor_z)), color = (:red, 0.3))
    
    # Create antenna (spherical cap above z=0)
    antenna_radius = ar.reactor_radius * 0.3
    antenna_center_z = ar.antenna_height / 2.0
    θ_antenna = collect(0:0.01:π)  # Upper semicircle
    antenna_x = antenna_radius .* cos.(θ_antenna)
    antenna_z = antenna_center_z .+ antenna_radius .* sin.(θ_antenna)
    lines!(ax, antenna_x, antenna_z, color = :blue, linewidth = 2, label = "Antenna")
    poly!(ax, Point2f.(zip(antenna_x, antenna_z)), color = (:blue, 0.3))
    
    # Mirror the antenna for x < 0
    lines!(ax, -antenna_x, antenna_z, color = :blue, linewidth = 2)
    
    ax.xlabel = "x (nm)"
    ax.ylabel = "z (nm)"
    ax.title = "Antenna-Reactor Structure"
    
    limit = ar.reactor_radius * 1.1
    xlims!(ax, -limit, limit)
    ylims!(ax, -ar.reactor_radius*1.1, ar.antenna_height*1.1)
    axislegend(ax, position = :rt)
    
    return fig
end

###Testing 
# Test CoreShell geometry
CS = build_CoreShell(25.0, 1.0, 5.0)
fig_cs = plot_coreshell_circles(CS)
display(fig_cs)

# Test CoreSatellite geometry
CSat = build_CoreSatellite(20.0, 5.0, 1.0)
fig_csat = plot_coresatellite_structure(CSat)
display(fig_csat)

# Test AntennaReactor geometry
AR = build_AntennaReactor(20.0, 30.0, 1.0)
fig_ar = plot_antennareactor_structure(AR)
display(fig_ar)
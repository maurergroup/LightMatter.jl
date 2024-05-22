using StaticArrays

abstract type Dimension end

struct Homogenous <: Dimension end

struct Linear{L} <: Dimension
    spatial_grid::SVector{L,Real}
end

struct Cylindrical{r,L} <: Dimension
    spatial_grid::SMatrix{r,L,Real}
end

struct Cubic{x,y,z} <: Dimension
    spatial_grid::SArray{Tuple{y,x,z},Real}
end

function define_sim_dimensions(sim::simulation_settings;Length::Integer,radius::Integer,
    xwidth::Integer,ywidth::Integer,spacing::Real)
    if sim.Dimension == 0
        return Homogenous
    elseif sim.Dimension == 1
        zs=0:spacing:Length
        return Linear{length(zs)}(zs)
    elseif sim.Dimension == 2
        rs=SVector{length(-radius:spacing:radius)}(-radius:spacing:radius)
        zs=SVector{length(0:spacing:Length)}(0:spacing:Length)
        slab_grid=Matrix{Tuple{Real,Real}}(undef,length(zs),length(rs))
        for (i,r) in enumerate(rs)
            for (j,z) in enumerate(zs)
                slab_grid[j,i]=(r,z)
            end
        end
        return Cylindrical{length(rs),length(zs)}(slab_grid)
    elseif sim.Dimension == 3
        xs=SVector{length(-xwidth/2:spacing:xwidth/2)}(-xwidth/2:spacing:xwidth/2)
        ys=SVector{length(-ywidth/2:spacing:ywidth/2)}(-ywidth/2:spacing:ywidth/2)
        zs=SVector{length(0:spacing:Length)}(0:spacing:Length)
        slab_grid=Array{Tuple{Real,Real,Real}}(undef,length(ys),length(xs),length(zs))
        for (i,x) in enumerate(xs)
            for (j,y) in enumerate(ys)
                for (k,z) in enumerate(zs)
                    slab_grid[j,i,k]=(x,y,z)
                end
            end
        end
        return Cubic{length(ys),length(xs),length(zs)}(slab_grid)
    end
end

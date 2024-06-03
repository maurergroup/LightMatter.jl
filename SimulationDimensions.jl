module SimulationDimensions
using StaticArrays
export Dimension, Homogenous, Linear, Cylindrical, Cubic, define_sim_dimensions
abstract type Dimension end

struct Homogenous <: Dimension end

@kwdef struct Linear{L} <: Dimension
    spatial_grid::SVector{L,Real}
end

@kwdef struct Cylindrical{r,L} <: Dimension
    spatial_grid::SMatrix{r,L,Real}
end

@kwdef struct Cubic{x,y,z} <: Dimension
    spatial_grid::SArray{Tuple{y,x,z},Real}
end

function define_sim_dimensions(;Dimension::Int64,Length=400::Integer,radius=200::Integer,
    xwidth=400::Integer,ywidth=400::Integer,spacing=1::Real)
    if Dimension == 0
        return Homogenous()
    elseif Dimension == 1
        return one_dimension(Length,spacing)
    elseif Dimension == 2
        return two_dimension(Length,spacing,radius)
    elseif Dimension == 3
        return three_dimension(length,xwidth,ywidth,spacing)
    end
end

function define_sim_dimensions(dict::Dict;Dimension=dict.Dimension::Int64,Length=dict.Length::Integer,radius=dict.Radius::Integer,
    xwidth=dict.XWidth::Integer,ywidth=dict.YWidth::Integer,spacing=dict.Spacing::Real)
    dict_overwrite(dict,Dimension=Dimension,Length=Length,radius=radius,xwidth=xwidth,ywidth=ywidth,spacing=spacing)
    if dict.Dimension == 0
        return Homogenous
    elseif dict.Dimension == 1
        return one_dimension(Length,spacing)
    elseif dict.Dimension == 2
        return two_dimension(Length,spacing,radius)
    elseif dict.Dimension == 3
        return three_dimension(length,xwidth,ywidth,spacing)
    end
end

function one_dimension(Length,spacing)
    zs=0:spacing:Length
    zs=collect(zs)
    slab = SVector{length(zs),typeof(zs[1])}(zs)
    return Linear{length(slab)}(spatial_grid=slab)
end

function two_dimension(Length,spacing,radius)
    rs=SVector{length(0:spacing:radius)}(0:spacing:radius)
    zs=SVector{length(0:spacing:Length)}(0:spacing:Length)
    slab_grid=Matrix{Tuple{Real,Real}}(undef,length(zs),length(rs))
    for (i,r) in enumerate(rs)
        for (j,z) in enumerate(zs)
            slab_grid[j,i]=(r,z)
        end
    end
    return Cylindrical{length(rs),length(zs)}(spatial_grid=slab_grid)
end

function three_dimension(Length,xwidth,ywidth,spacing)
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
    return Cubic{length(ys),length(xs),length(zs)}(spatial_grid=slab_grid)
end

end
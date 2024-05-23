using StaticArrays

abstract type Dimension end

struct Homogenous <: Dimension end

@kwdef struct Linear{L} <: Dimension
    Dim=1
    Boundaries=(L)
    spatial_grid::SVector{L,Real}
end

@kwdef struct Cylindrical{r,L} <: Dimension
    Dim=2
    Boundaries=(r,L)
    Spacing::Real
    spatial_grid::SMatrix{r,L,Real}
end

@kwdef struct Cubic{x,y,z} <: Dimension
    Dim=3
    Boundaries=(x,y,z)
    Spacing::Real
    spatial_grid::SArray{Tuple{y,x,z},Real}
end

function define_sim_dimensions(;Dimension::Int8,Length::Integer,radius::Integer,
    xwidth::Integer,ywidth::Integer,spacing::Real)
    if Dimension == 0
        return Homogenous
    elseif Dimension == 1
        return one_dimension(Length,spacing)
    elseif Dimension == 2
        return two_dimension(Length,spacing,radius)
    elseif Dimension == 3
        return three_dimension(length,xwidth,ywidth,spacing)
    end
end

function define_sim_dimensions(dict::Dict;Dimension=dict.Dimension::Int8,Length=dict.Length::Integer,radius=dict.Radius::Integer,
    xwidth=dict.XWidth::Integer,ywidth=dict.YWidth::Integer,spacing=dict.Spacing::Real)
    dict_overwrite(dict,Dimension=Dimension,Length=Length,radius=radius,xwidth=xwidth,ywidth=ywidth,spacing=spacing)
    if dict.Dimension == 0
        return Homogenous
    elseif dict.Dimension == 1
        return one_dimension(dict.Length,dict.spacing)
    elseif dict.Dimension == 2
        return two_dimension(dict.Length,dict.spacing,dict.radius)
    elseif dict.Dimension == 3
        return three_dimension(dict.length,dict.xwidth,dict.ywidth,dict.spacing)
    end
end

function one_dimension(length,spacing)
    zs=0:spacing:length
    return Linear{length(zs)}(Spacing=spacing,spatial_grid=zs)
end

function two_dimension(length,spacing,radius)
    rs=SVector{length(0:spacing:radius)}(0:spacing:radius)
    zs=SVector{length(0:spacing:Length)}(0:spacing:Length)
    slab_grid=Matrix{Tuple{Real,Real}}(undef,length(zs),length(rs))
    for (i,r) in enumerate(rs)
        for (j,z) in enumerate(zs)
            slab_grid[j,i]=(r,z)
        end
    end
    return Cylindrical{length(rs),length(zs)}(Spacing=spacing,spatial_grid=slab_grid)
end

function three_dimension(length,xwidth,ywidth,spacing)
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
    return Cubic{length(ys),length(xs),length(zs)}(Spacing=spacing,spatial_grid=slab_grid)
end

function dict_overwrite(dict::Dict;Dimension::Int8,Length::Integer,radius::Integer,
    xwidth::Integer,ywidth::Integer,spacing::Real)
    
end
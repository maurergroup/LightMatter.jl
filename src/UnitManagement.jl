#
#WIP!!! Currently not used anywhere
#
module UnitModule
using Unitful

@unit eVm "eVm" mass_for_eV (1/6.241509074460763e30)u"kg" true

end
"""
    convert_units(value::Union{Quantity, AbstractArray{<:Quantity}, Number, AbstractArray{<:Number}})

    Converts any user-given parameters that they have attached units to, to the correct units for Lightmatter.jl
    All values without Unitful.jl units attached are assumed to be in the correct units given by Lightmatter_units

    # Arguments
    - 'value': The value(array of values) of the parameter, either as a Unitful.jl Quantity or a Number

    # Returns
    - The Quantity converted into Lightmatter.jl's preferred units, or the Number left as is
"""
function convert_units(value::Union{Quantity, AbstractArray{<:Quantity}, Number, AbstractArray{<:Number}})
    if value isa Quantity || first(value) isa Quantity
        val = upreferred.(value)
        return Float64.(ustrip(val))
    else
        return value
    end
end

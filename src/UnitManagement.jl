#
#WIP!!! Currently not used anywhere
#
module UnitModule
using Unitful

@unit eVm "eVm" mass_for_eV (1/6.241509074460763e30)u"kg" true

end
"""
    convert_units(value::Union{Quantity, AbstractArray{<:Quantity}, Float64, AbstractArray{Float64}})

    Converts any user-given parameters that they have attached units to, to the correct units for LightMatter.jl
    All values without Unitful.jl units attached are assumed to be in the correct units given by LightMatter_units

    # Arguments
    - 'value': The value(array of values) of the parameter, either as a Unitful.jl Quantity or a Float64

    # Returns
    - The Quantity converted into LightMatter.jl's preferred units, or the Float64 left as is
"""
function convert_units(unit::Unitful.FreeUnits, value::Union{Quantity, AbstractArray{<:Quantity}, Float64, AbstractArray{Float64},Vector{Vector{Quantity}},Vector{Vector{Float64}}})
    if typeof(value) <: Vector{Vector{Quantity}}
        val = fill(zeros(length(value[1])), length(value))
        for i in eachindex(val)
            val[i] = Float64.(ustrip(uconvert.(unit, value[i])))
        end
        return val
    elseif value isa Quantity || first(value) isa Quantity
        val = uconvert.(unit, value)
        return Float64.(ustrip(val))
    else
        return value
    end
end

"""
    BaseUnits = 

    Global Tuple for conversion factors from base SI units to base LightMatter units
"""
global const BaseUnits = (time = 1e15, length = 1e9, mass = 6.2415e30, electric_current = 1, temperature = 1, amount = 1, luminosity = 1)

"""
    Constants = (ħ = 0.6582 eVfs, kB = 8.617e-5 eV/K, me = 5.686 eVm)

    Global named tuple for accessing constant physical values during a Simulation
"""
global const Constants = (ħ = ustrip(convert_units(u"eV*fs",Unitful.ħ)), kB = ustrip(convert_units(u"eV/K",Unitful.k)), me = ustrip(Unitful.me)*BaseUnits.mass, c = ustrip(convert_units(u"nm/fs",Unitful.c0)),
                          ϵ0 = 8.854e-12 * BaseUnits.time^4 / BaseUnits.mass/ BaseUnits.length^3, q=ustrip(-Unitful.q)*BaseUnits.time)


Unitful.uconvert(a::Unitful.FreeUnits, b::Union{Real, Array{<:Real}}) = b

"""
    LightMatter_units

    A list of the units used in LightMatter.jl: Please convert all units to this 
"""
global const LightMatter_units = [u"eV", u"nm", u"fs", u"K"] 


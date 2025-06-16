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
function convert_units(unit::Unitful.FreeUnits, value::Union{Quantity, AbstractArray{<:Quantity}, Number, AbstractArray{<:Number}})
    if value isa Quantity || first(value) isa Quantity
        val = uconvert.(unit, value)
        return Float64.(ustrip(val))
    else
        return value
    end
end

"""
    BaseUnits = 

    Global Tuple for conversion factors from base SI units to base Lightmatter units
"""
global const BaseUnits = (time = 1e15, length = 1e9, mass = 6.2415e30, electric_current = 1, temperature = 1, amount = 1, luminosity = 1)

"""
    Constants = (ħ = 0.6582 eVfs, kB = 8.617e-5 eV/K, me = 5.686 eVm)

    Global named tuple for accessing constant physical values during a Simulation
"""
global const Constants = (ħ = ustrip(convert_units(u"eV*fs",Unitful.ħ)), kB = ustrip(convert_units(u"eV/K",Unitful.k)), me = ustrip(Unitful.me)*BaseUnits.mass, c = ustrip(convert_units(u"nm/fs",Unitful.c0)),
                          ϵ0 = 8.854e-12 * BaseUnits.time^4 / BaseUnits.mass/ BaseUnits.length^3, q=ustrip(-Unitful.q)*BaseUnits.time)



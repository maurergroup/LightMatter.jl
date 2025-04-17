module UnitModule
using Unitful

@unit eVm "eVm" mass_for_eV (1/6.241509074460763e30)u"kg" true

end

function convert_units(value)
    if value isa Quantity || first(value) isa Quantity
        val = upreferred.(value)
        return Float64.(ustrip(val))
    else
        return value
    end
end

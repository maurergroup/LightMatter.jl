module UnitModule; using Unitful; @unit eVm "eVm" mass_for_eV (1/6.241509074460763e30)u"kg" true; end #Unit to make energy units eV as energy = m*L^2/t^2 and using L=fs and t = fs

function __init__()
    Unitful.register(UnitModule)
end

Unitful.register(UnitModule)

Unitful.preferunits(u"nm,fs,K,eVm"...)

function convert_units(value)
    if value isa Quantity || first(value) isa Quantity
        val = upreferred.(value)
        return ustrip(Float64.(val))
    else
        return value
    end
end

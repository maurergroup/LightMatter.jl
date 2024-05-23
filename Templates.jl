#=
    The role of this file is to set up some templates to populate the ODE's with for the various different ODE's. 
    This is necessary if we are to add BTE etc where there are no temperature ODE's present. Also ModelingToolkit
    cannot add terms to the right of an ODE just substitute temporary constants for variables
=#

abstract type Templates end

abstract type ElectronicTemperature <: Templates end
abstract type PhononicTemperature <: Templates end
abstract type ElectronDistribution <: Templates end

function ODEFactory()
    
end

function (::ElectronicTemperature)(;name)
    @parameters Spatial ElecPhon HeatCapacity
    @variables Tel(t) Source(t) 

    eqs = D(Tel) ~ (Source .+ Spatial .+ ElecPhon)./HeatCapacity

    ODESystem(eqs,t;name)
end

function (::PhononicTemperature)(;name)
    @parameters ElecPhon HeatCapacity
    @variables Tph(t) Source(t)  

    eqs = D(Tel) ~ (Source .+ ElecPhon)./HeatCapacity

    ODESystem(eqs,t;name)
end

function (::ElectronDistribution)(;name)
    @parameters Relax
    @variables Dis(t) Source(t)

    eqs = D(Dis) ~ Source .+ Relax

    ODESystem(eqs,t;name)

end
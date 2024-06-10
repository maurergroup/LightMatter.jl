function athem_builder(;named)
    @variables neqFD(t) Source(t) neqelel(t) neqelph(t)

    eqs = D(neqFD) ~ Source .+ neqelel .+ neqelph

    ODESystem(eqs,t;name)
end

function dFDdE(kB::Float64,Tel::Float64,μ::Float64,E::Float64)::Real
    Numer=-exp((E-μ)/(kB*Tel))
    Denom=kB*Tel*(exp((E-μ)/(kB*Tel))+1)^2
    return Numer/Denom
end

function dFDdT(kB::Float64,Tel::Float64,μ::Float64,E::Float64)::Real
    Numer=(E-μ)*exp((E-μ)/(kB*Tel))
    Denom=kB*Tel^2*(exp((E-μ)/(kB*Tel))+1)^2
    return Numer/Denom
end

FermiDirac(E::Float64,μ::Float64,Tel::Float64,kB::Float64) = 1/(exp((E-μ)/(kB*Tel))+1)

FermiDirac(E::Float64,μ::ForwardDiff.Dual,Tel::Float64,kB::Float64) = 1/(exp((E-μ)/(kB*Tel))+1)


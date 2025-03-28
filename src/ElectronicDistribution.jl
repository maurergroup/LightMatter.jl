"""
    FermiDirac(Tel::Real,μ::Real,E::Union{Vector{<:Real},Real})
    Returns the Fermi distribution at the given temperature and chemical potential for either a single energy
    value or across an energy window - given by whether E is a Real or Vector.
"""
@inline FermiDirac(Tel::Real,μ::Real,E::Union{Vector{<:Real},Real}) = 1 ./(exp.((E.-μ)./(Constants.kB*Tel)).+1)
"""
    dFDdE(Tel::Real,μ::Real,E::Union{Vector{<:Real},Real})
    Returns the derivative of a Fermi distribution with respect to energy at the current energy
    point or across an energy window - given by whether E is a Real or Vector. 
"""
@inline function dFDdE(Tel::Real,μ::Real,E::Union{Vector{<:Real},Real})
    return -exp.((E.-μ)./(Constants.kB*Tel))./(Constants.kB*Tel*(exp.((E.-μ)./(Constants.kB*Tel)).+1).^2)
end
"""
    dFDdT(Tel::Real,μ::Real,E::Union{Vector{<:Real},Real})
    Returns the derivative of a Fermi distribution with respect to temperature at the current energy
    point or across an energy window - given by whether E is a Real or Vector. 
"""
@inline function dFDdT(Tel::Real,μ::Real,E::Union{Vector{<:Real},Real})
    return (E.-μ).*exp.((E.-μ)./(Constants.kB*Tel))./(Constants.kB*Tel^2*(exp.((E.-μ)./(Constants.kB*Tel)).+1).^2)
end
"""
    dFDdμ(Tel::Real,μ::Real,E::Union{Vector{<:Real},Real})
    Returns the derivative of a Fermi distribution with respect to chemical potential at the current energy
    point or across an energy window - given by whether E is a Real or Vector. 
"""
@inline function dFDdμ(Tel::Real,μ::Real,E::Union{Vector{<:Real},Real})
    return exp.((E.-μ)./(Constants.kB*Tel))./(Constants.kB*Tel*(exp.((E.-μ)./(Constants.kB*Tel)).+1).^2)
end

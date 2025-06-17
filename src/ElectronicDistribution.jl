"""
    FermiDirac(Tel::Number, μ::Number, E::Union{Vector{<:Number},Number}) 
    
    Returns the thermal occupation of electrons at a given temperature, chemical potential and energy.
    If a vector of energies is given then it will return the distribution across that range.

    # Arguments
    - 'Tel': Thermal electronic temperature
    - 'μ': Chemical Potential
    - 'E': Energy value or range

    # Returns
    - Value or vector of the Fermi-Dirac distribution
"""
@inline FermiDirac(Tel, μ, E) = 1 ./ (exp.((E.-μ) ./ (Constants.kB*Tel)).+1)
"""
    dFDdE(Tel::Number, μ::Number, E::Union{Vector{<:Number},Number})
    
    Returns the change in the Fermi distribution with respect to energy at the given energy value or range.

    # Arguments
    - 'Tel': Thermal electronic temperature
    - 'μ': Chemical Potential
    - 'E': Energy value or range

    # Returns
    - Value or vector of the partial derivative of the Fermi distribution with respect to energy
"""
@inline function dFDdE(Tel, μ, E)
    return -exp.((E.-μ)./(Constants.kB*Tel)) ./ (Constants.kB*Tel * (exp.((E.-μ)./(Constants.kB*Tel)).+1).^2)
end
"""
    dFDdT(Tel::Number, μ::Number, E::Union{Vector{<:Number},Number})
    
    Returns the change in the Fermi distribution with respect to temperature at the given energy value or range.

    # Arguments
    - 'Tel': Thermal electronic temperature
    - 'μ': Chemical Potential
    - 'E': Energy value or range

    # Returns
    - Value or vector of the partial derivative of the Fermi distribution with respect to temperature
"""
@inline function dFDdT(Tel, μ, E)
    return (E.-μ) .* exp.((E.-μ)./(Constants.kB*Tel)) ./ (Constants.kB*Tel^2 * (exp.((E.-μ)./(Constants.kB*Tel)).+1).^2)
end
"""
    dFDdμ(Tel::Number, μ::Number, E::Union{Vector{<:Number},Number})
    
    Returns the change in the Fermi distribution with respect to chemical potential at the given energy value or range.

    # Arguments
    - 'Tel': Thermal electronic temperature
    - 'μ': Chemical Potential
    - 'E': Energy value or range

    # Returns
    - Value or vector of the partial derivative of the Fermi distribution with respect to chemical potential
"""
@inline function dFDdμ(Tel, μ, E)
    return exp.((E.-μ)./(Constants.kB*Tel)) ./ (Constants.kB*Tel * (exp.((E.-μ)./(Constants.kB*Tel)).+1).^2)
end

function magnetotransport_equations(sim)
    B = :($(sim.structure.fields.laser.magnetic) + $(sim.structure.fields.external.magnetic))
    return :(Lightmatter.magnetotransport_1d(fneq, sim, $B))
end

function df_dk(f, sim)
    k = sim.structure.bandstructure[2](sim.structure.egrid)
    fspl = DataInterpolations.AkimaInterpolation(f, k,extrapolation = ExtrapolationType.Constant)
    return DataInterpolations.derivative.(Ref(fspl), k)
end

function magnetotransport_1d(f, sim, B)

    h_2_e = findmin(abs.(sim.structure.egrid))[2]
    dfdk = Lightmatter.df_dk(f, sim)
    Y_e= Constants.q / Constants.ħ / Constants.c * sim.athermalelectrons.v_g[h_2_e+1:end] * B
    Y_h= -Constants.q / Constants.ħ / Constants.c * sim.athermalelectrons.v_g[1:h_2_e] * B

    tmp[1:h_2_e] = dfdk[1:h_2_e].*Y_h
    tmp[h_2_e+1:end] = dfdk[h_2_e+1:end].*Y_e
    return tmp * -1
end
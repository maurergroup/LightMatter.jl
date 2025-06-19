"""
    FermiDirac(Tel::Float64, μ::Float64, E::Union{Vector{Float64},Float64}) 
    
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
    dFDdE(Tel::Float64, μ::Float64, E::Union{Vector{Float64},Float64})
    
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
    dFDdT(Tel::Float64, μ::Float64, E::Union{Vector{Float64},Float64})
    
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
    dFDdμ(Tel::Float64, μ::Float64, E::Union{Vector{Float64},Float64})
    
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
    return :(Lightmatter.magnetotransport_1d!(Δf_mt, fneq, tot_n, sim, DOS, $B, g_k, tmp))#.+Lightmatter.FermiDirac(Tel, μ, sim.structure.egrid)
end

function df_dk!(dfdk::Vector{Float64}, g_k::Vector{Float64}, sim::Simulation)
    k = sim.structure.bandstructure[2](sim.structure.egrid)
    fspl = DataInterpolations.AkimaInterpolation(g_k, k, extrapolation = ExtrapolationType.Constant)

    @inbounds for i in eachindex(k)
        dfdk[i] = DataInterpolations.derivative(fspl, k[i])
    end

    return dfdk
end

function magnetotransport_1d!(Δf_mt::Vector{Float64}, f::Vector{Float64}, n::Float64, sim::Simulation, DOS::spl, B::Float64, g_k::Vector{Float64}, dfdk::Vector{Float64})
    h_2_e = get_h2e(sim)

    #= goal = Lightmatter.get_internalenergy(f, DOS, sim.structure.egrid)
    find_relaxeddistribution(g_k, sim.structure.egrid, goal, n, DOS)
    @inbounds @simd for i in eachindex(f)
        g_k[i] = f[i] - g_k[i]  # g_k now holds f - f₀
    end =#

    df_dk!(dfdk, f, sim)  # Compute derivative into dfdk (no aliasing)

    v_g = sim.athermalelectrons.v_g
    factor = Constants.q / (Constants.ħ * Constants.c) * B

    @inbounds @simd for i in 1:h_2_e
        Δf_mt[i] = -factor * v_g[i] * dfdk[i] * -1
    end
    @inbounds @simd for i in (h_2_e + 1):length(Δf_mt)
        Δf_mt[i] = factor * v_g[i] * dfdk[i] * -1
    end
end

function get_h2e(sim::Simulation)
    egrid = sim.structure.egrid
    h_2_e = 1
    best = abs(egrid[1])
    @inbounds for i in 2:length(egrid)
        v = abs(egrid[i])
        if v < best
            best = v
            h_2_e = i
        end
    end
    return h_2_e
end
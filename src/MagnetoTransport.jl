function magnetotransport_equations(sim)
    B = :($(sim.structure.fields.laser.magnetic) + $(sim.structure.fields.external.magnetic))
    return :(Lightmatter.magnetotransport_1d!(Δf_mt, fneq.+Lightmatter.FermiDirac(Tel, μ, sim.structure.egrid), sim, $B, DOS, n, band, sim.structure.egrid, g_k, tmp))
end

function df_dk!(dfdk::Vector{Float64}, f::Vector{Float64}, bandstructure::Vector{<:AkimaInterpolation}, egrid::Vector{Float64})
    k = bandstructure[2](egrid)
    fspl = DataInterpolations.AkimaInterpolation(f, k, extrapolation = ExtrapolationType.Constant)

    @inbounds for i in eachindex(k)
        dfdk[i] = DataInterpolations.derivative(fspl, k[i])
    end

    return dfdk
end

function magnetotransport_1d!(Δf_mt::Vector{Float64}, f::Vector{Float64}, sim::Simulation, B::Float64, DOS::spl, n::Float64, bandstructure::Vector{<:AkimaInterpolation}, egrid::Vector{Float64}, g_k::Vector{Float64}, dfdk::Vector{Float64})
    h_2_e = get_h2e(sim)

    goal = Lightmatter.get_internalenergy(f, DOS, sim.structure.egrid)
    find_relaxeddistribution(g_k, sim.structure.egrid, goal, n, DOS)
    @inbounds @simd for i in eachindex(f)
        g_k[i] = f[i] - g_k[i]  # g_k now holds f - f₀
    end

    df_dk!(dfdk, g_k, bandstructure, egrid)  # Compute derivative into dfdk (no aliasing)

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
function vonNeumann(dρ,ρ,H)
    dρ .= 1/(1im * Constants.ħ  )*(H*ρ - ρ*H)
    return nothing
end

function build_dm(N,occ)
    basis = []
    for i in 1:N
        ψ = zeros(N)
        ψ[i]=1
        push!(basis,ψ)
    end

    ρ_states = fill(zeros(N,N),N)
    for i in eachindex(basis)
        ρ_states[i] = occ[i]*basis[i]*basis[i]'
    end
    ρ = sum(ρ_states)
    return Complex.(ρ)
end

function thermal_bath_densitymatrix(H,β,ne)
    thermal_occupation = 1 ./(exp.(diag(H[1:end,1:end])*β).+1)
    ρ=Matrix(Diagonal(thermal_occupation))
    return Complex.(ρ./tr(ρ).*ne)
end

function discretize_DOS(dos_file, no_states, egrid)
    dos= readdlm(dos_file,comments=true)
    energy = dos[:,1]
    states = dos[:,2]
    DOS_tmp = generate_DOS(dos_file, 1.0)
    total_states = total_DOS_states(DOS_tmp, egrid)
    
    DOS = get_interpolant(energy, states*no_states ./ total_states)

    DOS_states = zeros(length(egrid))
    step = egrid[2]-egrid[1]
    for (i,E) in enumerate(egrid)
        int(u,p) = DOS(u)
        prob = IntegralProblem(int, E-step/2, E+step/2)
        DOS_states[i] = solve(prob, HCubatureJL(initdiv=1000), abstol=1e-8, reltol=1e-8).u
    end
    int_states  = round.(Int, DOS_states)
    disc_vector = vcat([fill(x,r) for (x,r) in zip(egrid, int_states)]...)
    return disc_vector
end

function total_DOS_states(DOS, egrid)
    int(u,p) = DOS(u)
    prob = IntegralProblem(int, egrid[1], egrid[end])
    return solve(prob, HCubatureJL(initdiv=1000), abstol=1e-8, reltol=1e-8).u
end

function construct_dipolevonNeumann(sim::Simulation)
    expr = quote
        println(t)
        sim = p.sim
        d̂ = sim.densitymatrix.DipoleMatrix
        H_tot = sim.densitymatrix.H0 .- (d̂[1] .* $(sim.densitymatrix.Fields.electric[1])) 
                                     .- (d̂[2] .* $(sim.densitymatrix.Fields.electric[2])) 
                                     .- (d̂[3] .* $(sim.densitymatrix.Fields.electric[3])) 
        Lightmatter.vonNeumann(du,u,H_tot)
    end
    return expr
end

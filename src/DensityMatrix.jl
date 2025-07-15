"""
    vonNeumann(dρ::Matrix{Float64},ρ::Matrix{Float64},H::Matrix{Float64})

    Equation of motion for a density matrix
    
    # Arguments
    - 'dρ': In-place derivative of the density matrix
    - 'ρ': Density matrix of the system
    - 'H': Hamiltonian matrix of the system

    # Returns
    - In-place function calculating dρ 
"""
function vonNeumann(dρ,ρ,H)
    dρ .= 1/(1im * Constants.ħ  )*(H*ρ - ρ*H)
    return nothing
end
"""
    build_dm(N::Int, occ::Vector{Float64})

    Builds a density matrix of NxN states with diagonal occupations occ
    
    # Arguments
    - 'N': Number of states
    - 'occ': Initial occupations of the density matrix 

    # Returns
    - An assembled density matrix with diagonal given by occ
"""
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
"""
    thermal_bath_densitymatrix(H::Matrix{Float64}, β::Float64, ne::Int)

    Builds a thermal density matrix assuming Fermi-Dirac statistics 
    
    # Arguments
    - 'H': Hamiltonian of the system, no external forces should be included here
    - 'Β': kB*T in the same units as the Hamiltonian
    - 'ne': The total number of electrons the density matrix represents

    # Returns
    - An assembled density matrix with a thermal occupation
"""
function thermal_bath_densitymatrix(H,β,ne)
    thermal_occupation = 1 ./(exp.(diag(H[1:end,1:end])*β).+1)
    ρ=Matrix(Diagonal(thermal_occupation))
    return Complex.(ρ./tr(ρ).*ne)
end
"""
    discretize_DOS(dos_file::String, no_states::Int, egrid::Vector{Float64})

    Discretizes a DOS and creates a vector of repeating energy states of length no_states.
    The number of these repeats should mimic the shape of the original DOS. The discretization is 
    only performed along the vector given by egrid. 
    
    # Arguments
    - 'dos_file': File path to DOS
    - 'no_states': Approximate length of the final vector (there may be slight differences in the final result)
    - 'egrid': Vector of energy values to discretize the DOS onto 

    # Returns
    - An approximation of a DOS to use as the diagonal of a Hamiltonian
"""
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
"""
    discretize_DOS(DOS::spl, egrid::Vector{Float64})

    Calculates the total number of states within the range of the egrid of the DOS
    provided by dos_file in discretize_DOS.
    
    # Arguments
    - 'DOS': Spline of the DOS generated from dos_file
    - 'egrid': Vector of energy values to discretize the DOS onto 

    # Returns
    - Total number of states in the energy window (egrid) of the DOS
"""
function total_DOS_states(DOS, egrid)::Float64
    int(u,p) = DOS(u)
    prob = IntegralProblem(int, egrid[1], egrid[end])
    return solve(prob, HCubatureJL(initdiv=1000), abstol=1e-8, reltol=1e-8).u
end
"""
    discretize_DOS(sim::Simulation)

    Constructs an expression for the propagation of a Hamiltonian's reaction to light
    in the dipole approximation. It takes the dipole matrix from sim.densitymatrix.DipoleMatrix
    
    # Arguments
    - 'sim': Settings of the simulation containing the dipole matrix and electric fields 

    # Returns
    - Expression for the propagation of a density matrix under illumination in the dipole approximation
"""
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

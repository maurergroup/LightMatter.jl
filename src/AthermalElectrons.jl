
"""
    athemdistribution_factory(sim, laser_expr)

    Constructs the time evolution of the non-equilibrium distribution in the AthEM model
    Find out more at https://arxiv.org/abs/2503.09479

    # Arguments
    - `sim`::Simulation
           Simulation struct containing physical and model parameters.
    - `laser_expr`::Expr
                  Expression representing the laser excitation term.

    # Returns
    - The total expression for the evolution of the athermal electron distribution, combining excitation and scattering terms.
"""
function athemdistribution_factory(sim::Simulation, laser_expr::Expr)
    feq = :(LightMatter.FermiDirac(Tel, μ, sim.structure.egrid))
    ftot = :($feq .+ fneq)
    Elecelec = athem_electronelectroninteraction(sim)
    Elecphon = athem_electronphononinteraction(sim)
    M = athemexcitation_matrixelements(sim)
    if sim.laser.hv isa Matrix
        athemexcite = quote
            Δfexcite .= vec(sum($laser_expr .* LightMatter.athemexcitation!(Δfexcite, tmp, $ftot, sim.structure.egrid, DOS, sim.laser.hv, $M, int_vec), dims=1))
            Δfexcite
        end
    else
        athemexcite = quote
                        LightMatter.athemexcitation!(Δfexcite, tmp, $ftot, sim.structure.egrid, DOS, sim.laser.hv, $M, int_vec, $laser_expr)
                        Δfexcite
                      end
    end
    return build_athemdistribution(sim, athemexcite, Elecelec, Elecphon)
end
"""
    build_athemdistribution(sim, athemexcite, Elecelec, Elecphon)

    Combines excitation and scattering contributions into a single broadcasted sum expression in the AthEM model

    # Arguments
    - `sim`::Simulation
           The simulation object.
    - `athemexcite`::Expr
                   Expression for athermal excitation.
    - `Elecelec`::Union{Expr, Float64}
                Electron-electron scattering Expr or 0.0 if disabled
    - `Elecphon`::Union{Expr, Float64}
                Electron-phonon scattering Expr or 0.0 if disabled

    # Returns
    - A broadcasted summation of all interaction terms.
"""
function build_athemdistribution(sim::Simulation, athemexcite::Expr, Elecelec::Union{Expr,Float64}, Elecphon::Union{Expr,Float64})
    args = Union{Expr, Symbol, Float64}[athemexcite, Elecelec, Elecphon]
    #= if sim.athermalelectrons.Conductivity == true
        push!(args, :f_cond)
    end =#
    return Expr(:call, :(.+), (args)...)
end
"""
    athemexcitation!(Δfneqe::DiffCache, Δfneqh::DiffCache, ftot::Vector{Float64}, egrid::Vector{Float64}, DOS::spl, hv::Float64, M::Union{Float64,Vector{Float64}})

    Computes the net non-equilibrium excitation using Fermi's Golden Rule.

    # Arguments
    - 'Δfneqe': A pre-initialized vector to reduce allocations. Can be DiffCache or Vector. 
                This vector contains the final result
    - 'Δfneqh': A pre-initialized vector to reduce allocations. Can be DiffCache or Vector.
    - `ftot`: Total distribution (f_eq + f_neq).
    - `egrid`: Energy grid.
    - `DOS`: Density of states spline.
    - `hv`: Photon energy.
    - `M`: Matrix elements.

    # Returns
    - In-place normalized change in particle distribution shape. Multiply by inputted laser energy at time t to get 
      full excitation shape change
"""
function athemexcitation!(Δfneqe, Δfneqh, ftot, egrid, DOS, hv::Float64, M, int_vec, laser)
    #Δfneqh = get_tmp(Δfneqh, ftot)
    #Δfneqe = get_tmp(Δfneqe, ftot)
    ftotspl = get_interpolant(egrid, ftot)
    athem_holegeneration!(Δfneqh, egrid, DOS, ftotspl, hv, M)
    athem_electrongeneration!(Δfneqe, egrid, DOS, ftotspl,hv, M)
    pc_sf = get_noparticles(int_vec, Δfneqe, DOS, egrid) / get_noparticles(int_vec, Δfneqh, DOS, egrid) # Corrects for particle conservation in the generation shape
    Δfneqe .-= (pc_sf * Δfneqh)

    #Δfneqe .*= excite_laser_internalenergy(Δfneqe, DOS, egrid, laser,int_vec)
    int_en = get_internalenergy(int_vec, Δfneqe, DOS, egrid)
    Δfneqe .*= ifelse(int_en == 0, laser, laser/int_en) # Normalizes the excitation to the inputted laser energy

    return nothing
end

function athemexcitation!(Δfneqh, Δfneqe, ftot, egrid, DOS, hv::Matrix{Float64}, M, int_vec,laser)
    Δneqs = zeros(size(hv,1),length(egrid))
    for i in 1:size(hv,1)
        ftotspl = get_interpolant(egrid, ftot)
        athem_holegeneration!(Δfneqh, egrid, DOS, ftotspl, hv[i,1], M)
        athem_electrongeneration!(Δfneqe, egrid, DOS, ftotspl,hv[i,1], M)
        pc_sf = get_noparticles(int_vec, Δfneqe, DOS, egrid) / get_noparticles(int_vec, Δfneqh, DOS, egrid) # Corrects for particle conservation in the generation shape
        Δfneqe .-= (pc_sf * Δfneqh)
        tot_en = get_internalenergy(int_vec, Δfneqe, DOS, egrid)
        Δneqs[i,:] .=  ifelse(tot_en == 0, Δfneqe .* hv[i,2], Δfneqe .* hv[i,2] ./ get_internalenergy(Δfneqe, DOS, egrid))
    end
    return Δneqs # Returns the different frequency changes
end

#= function excite_laser_internalenergy(dis, DOS, egrid, laser,int_vec)
    println(get_internalenergy(int_vec, 1.0*dis, DOS, egrid))
    int(u,p) = laser - get_internalenergy(int_vec, ForwardDiff.value(u)*dis, DOS, egrid)
    prob = NonlinearProblem(int, 1.0)
    sol = solve(prob, SimpleNewtonRaphson(), abstol=1e-3, reltol=1e-3).u
    println(sol)
    return sol
end =#
"""
    athem_holegeneration!(tmp::Vector, egrid::Vector{Float64},DOS::spl,ftotspl::spl,hv::Float64,M::Union{Float64,Vector{Float64}})

    Computes the shape of the distribution change due to hole generation.

    # Arguments
    - `tmp`: Temporary vector for in-place update
    - `egrid`: Energy grid.
    - `DOS`: Density of states spline.
    - `ftotspl::spl`: Spline of the total distribution.
    - `hv`: Photon energy.
    - `M`: Matrix elements.

    # Returns
    - Change in distribution due to hole generation
"""
function athem_holegeneration!(tmp, egrid, DOS, ftotspl, hv, M)
    for i in eachindex(egrid)
        tmp[i] = (2*pi/Constants.ħ) * M * DOS(egrid[i]+hv) * ftotspl(egrid[i]) * (1 - ftotspl(egrid[i]+hv))
    end
    return nothing
end
"""
    athem_electrongeneration!(tmp::Vector, egrid::Vector{Float64},DOS::spl,ftotspl::spl,hv::Float64,M::Union{Float64,Vector{Float64}})

    Computes the shape of the distribution change due to electron generation.

    # Arguments
    - `tmp`: Temporary vector for in-place update
    - `egrid`: Energy grid.
    - `DOS`: Density of states spline.
    - `ftotspl::spl`: Spline of the total distribution.
    - `hv`: Photon energy.
    - `M`: Matrix elements.

    # Returns
    - Change in distribution due to electron generation
"""
function athem_electrongeneration!(tmp, egrid, DOS, ftotspl, hv, M)
    for i in eachindex(egrid)
        tmp[i] = (2*pi/Constants.ħ) * M * DOS(egrid[i]-hv) * ftotspl(egrid[i]-hv) * (1 - ftotspl(egrid[i]))
    end
    return nothing
end
"""
    athemexcitation_matrixelements(sim::Simulation)

    Returns the symbolic matrix element expression for photo-excitation.
    Currently implemented:
    - :unity : No specified matrix elements - Returns 1.0

    # Arguments
    - 'sim': settings of the desired simulation

    # Returns
    - Matrix element expression. 
"""
function athemexcitation_matrixelements(sim::Simulation)
    if sim.athermalelectrons.ExcitationMatrixElements == :unity
        return :(1.0)
    end
end
"""
    athem_electronelectroninteraction(sim::Simulation)

    Returns an expression representing the electron-electron scattering term in the AthEM model.

    # Arguments
    - 'sim': settings of the desired simulation

    # Returns
    - Symbolic expression for relaxation or 0.0 if disabled.
"""
function athem_electronelectroninteraction(sim::Simulation)
    if sim.athermalelectrons.AthermalElectron_ElectronCoupling == true 
        expr = quote
            relax_dis
        end
        return expr # Uses relax_dis as a temp variable due to it being required here and in the electronic temperature system
    else
        return 0.0
    end
end
"""
    electron_relaxationtime(sim::Simulation)
    Returns the requested expression for the electronic relaxation time. 
    Currently implemented: 
    - :constant : A constant time defined by sim.athermalelectrons.τ
    - :FLT : Fermi Liquid lifetime defined in https://arxiv.org/abs/2503.09479

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Expression for the athermal electron lifetime due to electron-electron interactions.
"""
function electron_relaxationtime(sim::Simulation)
    if sim.athermalelectrons.ElectronicRelaxation == :constant
        return :(sim.athermalelectrons.τ)
    elseif sim.athermalelectrons.ElectronicRelaxation == :FLT
        return :(sim.athermalelectrons.τ * (μ-μ0+sim.athermalelectrons.FE)^2 ./((sim.structure.egrid.-μ).^2 .+ (pi*Constants.kB*Tel)^2))
    end
end       
"""
    athem_electronelectronscattering(fdis::VectorTel::Float64,μ::Float64,sim::Simulation,fneq::Vector{Float64},DOS::spl,n::Float64,τee::Union{Float64,Vector{Float64}})

    Calculates the electron-electron scattering contribution using a modified relaxation time approximation.

    # Arguments
    - 'Tel': Electronic temperature of the system
    - 'μ': Chemical potential at the current temperature
    - 'sim': Simulation settings and parameters
    - 'fneq': Current non-equilibrium electron distribution
    - 'DOS': Spline of the density-of-states
    - 'n': The number of electrons in the thermal system
    - 'τee': The athermal electron lifetime

    # Returns
    - Change in the non-equilibrium distribution due to scattering with a thermal electronic system
"""
function athem_electronelectronscattering!(fdis, frel, Tel, μ, egrid, fneq, DOS, n, τee, int_vec, tmp)
    ftot = LightMatter.FermiDirac(Tel, μ, egrid) .+ fneq
    goal = get_internalenergy(int_vec, ftot, DOS, egrid)
    find_relaxeddistribution(frel, egrid, goal, n, DOS, int_vec, μ, tmp)
    fdis .= -(fneq ./τee) .+ ((ftot-frel) ./ τee) #fdis = feq + fneq so we need to add fneq twice to get fneq - feq
end
"""
    find_relaxeddistribution(egrid::Vector{Float64},goal::Float64,n::Float64,DOS::spl,kB::Float64)
    
    Solves for a Fermi distribution with the same internal energy as a given target.

    # Arguments
    - 'egrid': Energy grid distributions are solved on
    - 'goal': The internal energy of the total electronic system
    - 'n': Float64 of electrons in the thermal system
    - 'DOS': Spline of the density-of-states

    # Returns
    - Fermi-Dirac distribution with same internal energy as the goal.
"""
function find_relaxeddistribution(out, egrid, goal, n, DOS, int_vec, μ0, tmp)
    prob = IntervalNonlinearProblem(find_relaxedtemp, (100.0, 1e5), (out, n, DOS, egrid, goal, int_vec, tmp, μ0))
    sol = solve(prob; alg=Brent(),abstol=1e-5, reltol=1e-5).u
    μ = find_chemicalpotential(n, sol, DOS, egrid, int_vec, tmp, μ0)
    #out = get_tmp(out, sol)
    FermiDirac!(out, sol, μ, egrid)
    return nothing
end
"""
    find_relaxedtemp(Tel::Float64,n::Float64,DOS::spl,egrid::Vector{Float64})

    Given a temperature guess, computes chemical potential and internal energy.

    # Arguments
    - 'Tel': Guessed electronic temperature to match the goal
    - 'n': Float64 of electrons in the thermal system
    - 'egrid': Energy grid distributions are solved on
    - 'DOS': Spline of the density-of-states

    # Returns
    - Internal energy of the current temperature guess.
"""
function find_relaxedtemp(u, (out, n, DOS, egrid, goal, int_vec, tmp, μ0))
    μ = find_chemicalpotential(n, u, DOS, egrid, int_vec, tmp, μ0)
    FermiDirac!(out, u, μ, egrid)
    return goal - get_internalenergy(int_vec, out, DOS, egrid)
end 
"""
    athem_electronelectroninteraction(sim::Simulation)

    Returns an expression representing the athermal electron-phonon scattering term in the AthEM model.

    # Arguments
    - 'sim': settings of the desired simulation

    # Returns
    - Symbolic expression for relaxation or 0.0 if disabled.
"""
function athem_electronphononinteraction(sim::Simulation)
    if sim.athermalelectrons.AthermalElectron_PhononCoupling == true 
        τep = phonon_relaxationtime(sim)
        return :(-fneq / $τep)
    else
        return 0.0
    end
end
"""
    phonon_relaxationtime(sim::Simulation)

    Determines the requested expression for the electron-phonon relaxation time. 
    Currently implemented:
    - :constant : a constant lifetime given by simulation.athermalelectrons.τep
    - :quasi : quasiparticle relaxation time defined in E.Carpene, Phys. Rev. B., 2006, 74, 024301.

    # Arguments
    - 'sim': settings of the desired simulation

    # Returns
    - Expression for the athermal electron lifetime due to phonon interactions.
"""
function phonon_relaxationtime(sim::Simulation)
    if sim.athermalelectrons.PhononicRelaxation == :constant
        return :(sim.athermalelectrons.τep)
    elseif sim.athermalelectrons.PhononicRelaxation == :quasi
        return :(sim.electronictemperature.γ * (Tel+Tph) / (2*sim.electronictemperature.g))
    end
end  
"""
    athem_thermalelectronparticlechange(sim::Simulation)
    
    Returns an expression for the total thermal electronum number due to relaxation and optional conductivity.

    # Arguments
    - 'sim': settings of the desired simulation

    # Returns
    - Expr for the time dependence of the thermal electron number.
"""
function athem_thermalelectronparticlechange(sim::Simulation)
    return :(-LightMatter.get_noparticles(int_vec, du.fneq[i,:],DOS,sim.structure.egrid))
end
"""
    electron_distribution_transport!(v_g::Vector{Float64},f::AbstractArray{Float64},Δf::AbstractArray{Float64},dz::Float64)

    Computes ballistic transport of an electronic distribution using second-order finite differences via a kinetic like model.
    Uses forward(reverse) difference for the boundaries. 

    # Arguments
    - `v_g`: Group velocity vector or vector of vectors(spatially-resolved DOS)
    - `f`: Distribution function.
    - `Δf`: Output array to store result.
    - `dz`: Spatial resolution.

    # Returns
    - In-place change to Δf
"""
function electron_distribution_transport!(Δf, v_g::Vector{Float64}, f, dz, Tel, noe, ftot, tmp, sim)
    calculate_ftot(f, Tel, noe, ftot, tmp, sim)
    
    @views @inbounds for i in 2:size(f, 1)-1
        # Ballistic transport: advection in space
        @. Δf[i,:] = -v_g * (ftot[i+1,:] - ftot[i-1,:]) / (2*dz)
    end

    @views @. Δf[1, :] = -v_g * (ftot[2, :] - ftot[1, :]) / dz
    @views @. Δf[end, :] = -v_g * (ftot[end, :] - ftot[end-1, :]) / dz
end

function electron_distribution_transport!(Δf, v_g::Matrix{Float64}, f, dz, Tel, noe, ftot, tmp, sim)
    calculate_ftot(f, Tel, noe, ftot, tmp, sim)
    
    @views @inbounds for i in 2:size(f, 1)-1
      # Ballistic transport: advection in space
      @. Δf[i,:] = -v_g[i,:] * (ftot[i+1,:] - ftot[i-1,:]) / (2*dz)
    end
    @views @. Δf[1,:] = -v_g[1,:] * (ftot[2,:] - ftot[1,:]) / dz
    @views @. Δf[end,:] = -v_g[end,:] * (ftot[end,:] - ftot[end-1,:]) / dz
end

function calculate_ftot(f, Tel::Real, noe, tmp, tmp2, sim)
    if sim.structure.Elemental_System > 1
        @views @inbounds for i in 1:size(f, 1)
            X = LightMatter.mat_picker(sim.structure.dimension.grid[i], sim.structure.dimension.InterfaceHeight)
            μ = LightMatter.find_chemicalpotential(noe[i], Tel, sim.structure.DOS[X], sim.structure.egrid, tmp[i,:], tmp2[i,:],sim.structure.μ_offset[X])
            LightMatter.FermiDirac!(view(tmp, i, :),Tel, μ, sim.structure.egrid)
            tmp[i,:] .+= f[i,:]
        end
    else
        @views @inbounds for i in 1:size(f, 1)
            μ = LightMatter.find_chemicalpotential(noe[i], Tel, sim.structure.DOS, sim.structure.egrid, tmp[i,:], tmp2[i,:], sim.structure.μ_offset[1])
            LightMatter.FermiDirac!(view(tmp, i, :),Tel, μ, sim.structure.egrid)
            tmp[i,:] .+= f[i,:]
        end
    end
end

function calculate_ftot(f, Tel::AbstractVector, noe, tmp, tmp2, sim)
    if sim.structure.Elemental_System > 1
        @views @inbounds for i in 1:size(f, 1)
            X = LightMatter.mat_picker(sim.structure.dimension.grid[i], sim.structure.dimension.InterfaceHeight)
            μ = LightMatter.find_chemicalpotential(noe[i], Tel[i], sim.structure.DOS[X], sim.structure.egrid, tmp[i,:], tmp2[i,:], sim.structure.μ_offset[X])
            LightMatter.FermiDirac!(view(tmp, i, :),Tel[i], μ, sim.structure.egrid)
            @views tmp[i,:] .+= f[i,:]
        end
    else
        @views @inbounds for i in 1:size(f, 1)
            μ = LightMatter.find_chemicalpotential(noe[i], Tel[i], sim.structure.DOS, sim.structure.egrid, tmp[i,:], tmp2[i,:], sim.structure.μ_offset[1])
            LightMatter.FermiDirac!(view(tmp, i, :),Tel[i], μ, sim.structure.egrid)
            tmp[i,:] .+= f[i,:]
        end
    end
end

"""
    FE_initalization(bulk_DOS::Union{String, Vector{String}})

    Initializes Fermi energy based on bulk density of states. This is defined as the difference between the bottom of the
    valence band and 0.0. This is required for calculations such as the ballistic velocity or athermal electron lifetime due
    to electronic interactions.

    # Arguments
    - `bulk_DOS`: (Array of) file path to bulk density-of-states

    # Returns
    - Fermi energy/energies for each DOS input.
"""
function FE_initialization(bulk_DOS::String)
    FE = get_FermiEnergy(bulk_DOS)
    return FE
end

function create_energy_window(egrid::Vector{Float64}, center::Float64; width_below::Float64=0.05, width_above::Float64=0.05, fill_value::Float64=1.0)
    result = zeros(Float64, length(egrid))
    lower_bound = center - width_below
    upper_bound = center + width_above
    
    @inbounds for i in eachindex(egrid)
        if lower_bound <= egrid[i] <= upper_bound
            result[i] = fill_value
        end
    end
    
    return result
end

function FE_initialization(bulk_DOS::Vector{String}, μ_offset::Bool=true, reference=1)
    FE = zeros(length(bulk_DOS))
    if μ_offset
        offset = calculate_μoffset(bulk_DOS, reference)
        for i in eachindex(bulk_DOS)
            FE[i] = get_FermiEnergy(bulk_DOS[i]) + offset[i]
        end
    else
        for i in eachindex(bulk_DOS)
            FE[i] = get_FermiEnergy(bulk_DOS[i])
        end
    end
    return FE
end

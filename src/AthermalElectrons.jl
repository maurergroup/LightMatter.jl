
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
            Δfexcite .= vec(sum($laser_expr .* LightMatter.athemexcitation!(Δfexcite, tmp, $ftot, sim.structure.egrid, DOS, sim.laser.hv, $M), dims=1))
            Δfexcite
        end
    else
        athemexcite = quote
                        LightMatter.athemexcitation!(Δfexcite, tmp, $ftot, sim.structure.egrid, DOS, sim.laser.hv, $M)
                        Δfexcite = LightMatter.access_DiffCache(Δfexcite, fneq[1])
                        $laser_expr .* Δfexcite
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
    if sim.athermalelectrons.Conductivity == true
        push!(args, :f_cond)
    end
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
function athemexcitation!(Δfneqe, Δfneqh, ftot, egrid, DOS, hv::Float64, M)
    Δfneqh = get_tmp(Δfneqh, ftot)
    Δfneqe = get_tmp(Δfneqe, ftot)
    ftotspl = get_interpolant(egrid, ftot)
    athem_holegeneration!(Δfneqh, egrid, DOS, ftotspl, hv, M)
    athem_electrongeneration!(Δfneqe, egrid, DOS, ftotspl,hv, M)
    pc_sf = get_noparticles(Δfneqe, DOS, egrid) / get_noparticles(Δfneqh, DOS, egrid) # Corrects for particle conservation in the generation shape
    Δfneqe .-= (pc_sf * Δfneqh)
    Δfneqe ./= get_internalenergy(Δfneqe, DOS, egrid) # Scales the shape of the change by the internal energy to later match with the laser
    return nothing
end

function athemexcitation!(Δfneqh, Δfneqe, ftot, egrid, DOS, hv::Matrix{Float64}, M)
    Δneqs = zeros(size(hv,1),length(egrid))
    for i in 1:size(hv,1)
        ftotspl = get_interpolant(egrid, ftot)
        athem_holegeneration!(Δfneqh, egrid, DOS, ftotspl, hv[i,1], M)
        athem_electrongeneration!(Δfneqe, egrid, DOS, ftotspl,hv[i,1], M)
        pc_sf = get_noparticles(Δfneqe, DOS, egrid) / get_noparticles(Δfneqh, DOS, egrid) # Corrects for particle conservation in the generation shape
        Δfneqe .-= (pc_sf * Δfneqh)
        Δneqs[i,:] .= Δfneqe .* hv[i,2] ./ get_internalenergy(Δfneqe, DOS, egrid) # Scales by internal energy and fraction of fluence at given frequency (hv[i][2])
    end
    return Δneqs # Returns the different frequency changes
end
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
            relax_dis = LightMatter.access_DiffCache(relax_dis, fneq[1])
            -1 * relax_dis
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
        return :(sim.athermalelectrons.τ * (μ+sim.athermalelectrons.FE)^2 ./((sim.structure.egrid.-μ).^2 .+ (pi*Constants.kB*Tel)^2))
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
function athem_electronelectronscattering!(fdis, frel, Tel, μ, egrid, fneq, DOS, n, τee)
    fdis = get_tmp(fdis, Tel)
    LightMatter.FermiDirac!(fdis, Tel, μ, egrid)
    fdis .+= fneq
    goal = get_internalenergy(fdis, DOS, egrid)
    find_relaxeddistribution!(frel, egrid, goal, n, DOS)
    frel = get_tmp(frel, Tel)
    @. fdis = (fneq+frel-fdis+fneq) ./ τee #fdis = feq + fneq so we need to subtract fneq twice to get feq - fneq
    return nothing
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
function find_relaxeddistribution!(out, egrid, goal, n, DOS)
    prob = NonlinearProblem(find_temperatureandμ!, SA[1000.0,0.0], (out, n, DOS, egrid, goal))
    sol = solve(prob, SimpleNewtonRaphson(); abstol=1e-10, reltol=1e-10).u
    out = get_tmp(out, n)
    FermiDirac!(out, sol[1], sol[2], egrid)
    return nothing
end
"""
    find_temperatureandμ(Tel::Float64,n::Float64,DOS::spl,egrid::Vector{Float64})

    Given a temperature guess, computes chemical potential and internal energy.

    # Arguments
    - 'Tel': Guessed electronic temperature to match the goal
    - 'n': Float64 of electrons in the thermal system
    - 'egrid': Energy grid distributions are solved on
    - 'DOS': Spline of the density-of-states

    # Returns
    - Internal energy of the current temperature guess.
"""
function find_temperatureandμ!(du, u, (out, n, DOS, egrid, goal))
    n = ForwardDiff.value(n)
    goal = ForwardDiff.value(goal)
    out = get_tmp(out, u[1])
    FermiDirac!(out,u[1], u[2],egrid)
    du[1] = goal - get_internalenergy(out, DOS, egrid)
    du[2] = n - get_noparticles(out, DOS, egrid)
    return nothing
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
    return :(-LightMatter.get_noparticles(du.fneq[i,:],DOS,sim.structure.egrid))
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
function electron_distribution_transport!(v_g::Vector{Float64}, f::Matrix{Float64}, Δf::Matrix{Float64}, dz::Real)
    @views @inbounds for i in 2:size(f, 1)-1
        @. Δf[i, :] = ((f[i-1, :] - 2 * f[i, :] + f[i+1, :]) / dz) * v_g
    end

    @views @. Δf[1, :] = (-(f[1, :] - f[2, :]) / dz) * v_g
    @views @. Δf[end, :] = ((f[end-1, :] - f[end, :]) / dz) * v_g
end

function electron_distribution_transport!(v_g::Matrix{Float64}, f, Δf, dz)
    @views @inbounds for i in 2:size(f, 1)-1
        @. Δf[i,:] = (f[i-1,:] - 2*f[i,:] + f[i+1,:]) / dz * v_g[i,:]
    end
    @views @. Δf[1,:] = -(f[1,:] - f[2,:]) ./ dz *v_g[1,:]
    @views @. Δf[end,:] = (f[end-1,:] - f[end,:]) / dz * v_g[end,:]
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
function FE_initialization(bulk_DOS::Union{String, Vector{String}})
    if bulk_DOS isa String
        FE = get_FermiEnergy(bulk_DOS)
        return FE
    else
        FE = zeros(length(bulk_DOS))
        for i in eachindex(bulk_DOS)
            FE[i] = get_FermiEnergy(bulk_DOS[i])
        end
        return FE
    end
end

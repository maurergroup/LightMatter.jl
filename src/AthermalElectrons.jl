"""
    athemdistribution_factory(sim::Simulation, laser_expr::Expr)

    Constructs the time evolution of the non-equilibrium distribution in the AthEM model
    Find out more at https://arxiv.org/abs/2503.09479

    # Arguments
    - `sim`: Simulation struct containing physical and model parameters.
    - `laser_expr`: Expression representing the laser excitation term.

    # Returns
    - The total expression for the evolution of the athermal electron distribution, combining excitation and scattering terms.
"""
function athemdistribution_factory(sim::Simulation, laser_expr::Expr)
    feq = :(Lightmatter.FermiDirac(Tel, μ, sim.structure.egrid))
    ftot = :($feq .+ fneq)
    Elecelec = athem_electronelectroninteraction(sim)
    Elecphon = athem_electronphononinteraction(sim)
    mag_trans = athem_magneotransport(sim)
    M = athemexcitation_matrixelements(sim)
    if sim.laser.hv isa Matrix
        athemexcite = :(vec(sum($laser_expr .* Lightmatter.athemexcitation($ftot, sim.structure.egrid, DOS, sim.laser.hv, $M), dims=1)))
    else
        athemexcite = :($laser_expr .* Lightmatter.athemexcitation($ftot, sim.structure.egrid, DOS, sim.laser.hv, $M))
    end
    return build_athemdistribution(sim, athemexcite, Elecelec, Elecphon, mag_trans)
end
"""
    build_athemdistribution(sim::Simulation, athemexcite::Expr, Elecelec::Union{Expr, Number}, Elecphon::Union{Expr, Number})

    Combines excitation and scattering contributions into a single broadcasted sum expression in the AthEM model

    # Arguments
    - `sim`: The simulation object.
    - `athemexcite`: Expression for athermal excitation.
    - `Elecelec`: Electron-electron scattering Expr or 0.0 if disabled
    - `Elecphon`: Electron-phonon scattering Expr or 0.0 if disabled

    # Returns
    - A broadcasted summation of all interaction terms.
"""
function build_athemdistribution(sim::Simulation, athemexcite::Expr, Elecelec::Union{Expr,Number}, Elecphon::Union{Expr,Number}, mag_trans::Union{Symbol,Number})
    args = Union{Expr, Symbol, Number}[athemexcite, Elecelec, Elecphon, mag_trans]
    if sim.athermalelectrons.Conductivity == true
        push!(args, :f_cond)
    end
    return Expr(:call, :(.+), (args)...)
end
"""
    athemexcitation(ftot::Vector{<:Number}, egrid::Vector{<:Number}, DOS::spl, hv::Number, M::Union{Number,Vector{<:Number}})

    Computes the net non-equilibrium excitation using Fermi's Golden Rule.

    # Arguments
    - `ftot`: Total distribution (f_eq + f_neq).
    - `egrid`: Energy grid.
    - `DOS`: Density of states spline.
    - `hv`: Photon energy.
    - `M`: Matrix elements.

    # Returns
    - Normalized excitation change in the distribution.
"""
function athemexcitation(ftot, egrid, DOS, hv::Number, M)
    ftotspl = get_interpolant(egrid, ftot)
    Δfneqh = athem_holegeneration(egrid, DOS, ftotspl, hv, M)
    Δfneqe = athem_electrongeneration(egrid, DOS, ftotspl,hv, M)
    pc_sf = get_noparticles(Δfneqe, DOS, egrid) / get_noparticles(Δfneqh, DOS, egrid) # Corrects for particle conservation in the generation shape
    Δfneqtot = Δfneqe .- (pc_sf * Δfneqh)
    return Δfneqtot ./ get_internalenergy(Δfneqtot, DOS, egrid) # Scales the shape of the change by the internal energy to later match with the laser
end

function athemexcitation(ftot, egrid, DOS, hv::Matrix{<:Number}, M)
    Δneqs = zeros(size(hv,1),length(egrid))
    for i in 1:size(hv,1)
        ftotspl = get_interpolant(egrid, ftot)
        Δfneqh = athem_holegeneration(egrid, DOS, ftotspl, hv[i,1], M)
        Δfneqe = athem_electrongeneration(egrid, DOS, ftotspl,hv[i,1], M)
        pc_sf = get_noparticles(Δfneqe, DOS, egrid) / get_noparticles(Δfneqh, DOS, egrid) # Corrects for particle conservation in the generation shape
        Δfneqtot = Δfneqe .- (pc_sf * Δfneqh)
        Δneqs[i,:] .= Δfneqtot .* hv[i,2] ./ get_internalenergy(Δfneqtot, DOS, egrid) # Scales by internal energy and fraction of fluence at given frequency (hv[i][2])
    end
    return Δneqs # Returns the different frequency changes
end
"""
    athem_holegeneration(egrid::Vector{<:Number},DOS::spl,ftotspl::spl,hv::Number,M::Union{Number,Vector{<:Number}})

    Computes the shape of the distribution change due to hole generation.

    # Arguments
    - `egrid`: Energy grid.
    - `DOS`: Density of states spline.
    - `ftotspl::spl`: Spline of the total distribution.
    - `hv`: Photon energy.
    - `M`: Matrix elements.

    # Returns
    - Change in distribution due to hole generation
"""
function athem_holegeneration(egrid, DOS, ftotspl, hv, M)
    return (2*pi/Constants.ħ) .* M .* DOS(egrid.+hv) .* ftotspl(egrid) .* (1 .- ftotspl(egrid.+hv))
end
"""
    athem_electrongeneration(egrid::Vector{<:Number},DOS::spl,ftotspl::spl,hv::Number,M::Union{Number,Vector{<:Number}})

    Computes the shape of the distribution change due to electron generation.

    # Arguments
    - `egrid`: Energy grid.
    - `DOS`: Density of states spline.
    - `ftotspl::spl`: Spline of the total distribution.
    - `hv`: Photon energy.
    - `M`: Matrix elements.

    # Returns
    - Change in distribution due to electron generation
"""
function athem_electrongeneration(egrid, DOS, ftotspl, hv, M)
    return (2*pi/Constants.ħ) .* M .* DOS(egrid.-hv) .* ftotspl(egrid.-hv) .* (1 .-ftotspl(egrid))
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
        return :(-1 * relax_dis) # Uses relax_dis as a temp variable due to it being required here and in the electronic temperature system
    else
        return 0.0
    end
end
"""
    athem_electronelectronscattering(Tel::Number,μ::Number,sim::Simulation,fneq::Vector{<:Number},DOS::spl,n::Number,τee::Union{Number,Vector{<:Number}})

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
function athem_electronelectronscattering(Tel, μ, sim::Simulation, fneq, DOS, n, τee)
    feq = Lightmatter.FermiDirac(Tel, μ, sim.structure.egrid)
    ftot = feq .+ fneq
    goal = extended_Bode(ftot.*DOS(sim.structure.egrid).*sim.structure.egrid, sim.structure.egrid)
    frel = find_relaxeddistribution(sim.structure.egrid, goal, n, DOS)
    return (fneq.+frel.-feq) ./ τee
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
        return :(sim.athermalelectrons.τ * (μ.+sim.athermalelectrons.FE)^2 ./((sim.structure.egrid.-μ).^2 .+ (pi*Constants.kB*Tel)^2))
    end
end       
"""
    find_relaxeddistribution(egrid::Vector{<:Number},goal::Number,n::Number,DOS::spl,kB::Number)
    
    Solves for a Fermi distribution with the same internal energy as a given target.

    # Arguments
    - 'egrid': Energy grid distributions are solved on
    - 'goal': The internal energy of the total electronic system
    - 'n': Number of electrons in the thermal system
    - 'DOS': Spline of the density-of-states

    # Returns
    - Fermi-Dirac distribution with same internal energy as the goal.
"""
function find_relaxeddistribution(egrid, goal, n, DOS)
    f(u,p) = goal - find_temperatureandμ(u, n, DOS, egrid)
    Temp = solve(NonlinearProblem(f,1000.0); abstol=1e-3, reltol=1e-3).u
    μ = find_chemicalpotential(n, Temp, DOS, egrid)
    return FermiDirac(Temp, μ, egrid)
end

#= function find_relaxeddistribution(egrid, goal::ForwardDiff.Dual, noe::ForwardDiff.Dual, DOS)
    int = ForwardDiff.value(goal)
    n = ForwardDiff.value(noe)
    f(u,p) = int - find_temperatureandμ(u, n, DOS, egrid)
    Temp = solve(NonlinearProblem(f,1000.0); abstol=1e-10, reltol=1e-10).u
    μ = find_chemicalpotential(n, Temp, DOS, egrid)
    return FermiDirac(Temp, μ, egrid)
end =#
"""
    find_temperatureandμ(Tel::Number,n::Number,DOS::spl,egrid::Vector{<:Number})

    Given a temperature guess, computes chemical potential and internal energy.

    # Arguments
    - 'Tel': Guessed electronic temperature to match the goal
    - 'n': Number of electrons in the thermal system
    - 'egrid': Energy grid distributions are solved on
    - 'DOS': Spline of the density-of-states

    # Returns
    - Internal energy of the current temperature guess.
"""
function find_temperatureandμ(Tel, n, DOS, egrid)
    μ = find_chemicalpotential(n, Tel, DOS, egrid)
    return get_internalenergy(FermiDirac(Tel,μ,egrid), DOS, egrid)
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
        return :(-fneq ./ $τep)
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
    part_change = :(Lightmatter.get_noparticles(relax_dis, DOS, sim.structure.egrid))
    args = Vector{Union{Symbol,Expr}}([part_change])
    if sim.athermalelectrons.Conductivity == true
        push!(args, :n_cond)
    end
    if sim.athermalelectrons.AthermalElectron_PhononCoupling == true
        τep = phonon_relaxationtime(sim)
        push!(args, :(Lightmatter.get_noparticles(fneq./$τep,DOS,egrid)))
    end
    if sim.athermalelectrons.MagnetoTransport == true
        push!(args, :(Lightmatter.get_noparticles(Δf_mt,DOS,egrid)))
    end
    return Expr(:call, :+, args...)
end
"""
    electron_distribution_transport!(v_g::Vector{<:Number},f::AbstractArray{<:Number},Δf::AbstractArray{<:Number},dz::Number)

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
function electron_distribution_transport!(v_g::Vector{<:Number}, f, Δf, dz)
    for i in 2:size(f, 1)-1
        Δf[i,:] = (f[i-1,:] .- 2*f[i,:] .+ f[i+1,:]) ./ dz .* v_g
    end
    Δf[1,:] = -(f[1,:] .- f[2,:]) ./ dz .*v_g
    Δf[end,:] = (f[end-1,:] .- f[end,:]) ./ dz .* v_g
end

function electron_distribution_transport!(v_g::Matrix{<:Number}, f, Δf, dz)
    for i in 2:size(f, 1)-1
        Δf[i,:] = (f[i-1,:] .- 2*f[i,:] .+ f[i+1,:]) ./ dz .* v_g[i,:]
    end
    Δf[1,:] = -(f[1,:] .- f[2,:]) ./ dz .*v_g[1,:]
    Δf[end,:] = (f[end-1,:] .- f[end,:]) ./ dz .* v_g[end,:]
end
"""
    thermal_particle_transport(v_g::Vector{<:Number},egrid::Vector{<:Number},n::Vector{<:Number},Δn::Vector{<:Number},dz::Number)

    Calculates transport correction for thermal particle distribution using the Fermi velocity.

    # Arguments
    - `v_g`: Group velocity vector or vector of vectors(spatially-resolved DOS)
    - `egrid`: Energy grid.
    - `n`: Number of particles in thermal system
    - `Δn`: Output array.
    - `dz`: Spatial resolution.

    # Returns
    - Updated transport correction array.
"""
function thermal_particle_transport!(v_g::Vector{<:Number}, egrid, n, Δn, dz)
    idx_0 = findmin(abs.(egrid.-0.0))[2]
    v_F = v_g[idx_0]
    for i in 2:length(Δn) -1
        Δn[i] = (n[i+1] - (2*n[i]) + n[i-1]) / dz * v_F
    end
    Δn[1] = (n[2] - n[1]) / dz * v_F
    Δn[end] = (n[end-1] - n[end]) / dz * v_F
end

function thermal_particle_transport!(v_g::Matrix{<:Number}, egrid, n, Δn, dz)
    idx_0 = findmin(abs.(egrid.-0.0))[2]
    for i in 2:length(Δn) -1
        v_F = v_g[i,idx_0]
        Δn[i] = (n[i+1] - (2*n[i]) + n[i-1]) / dz * v_F
    end
    v_F1 = v_g[1,idx_0]
    Δn[1] = (n[2] - n[1]) / dz * v_F1
    v_Fend = v_g[end,idx_0]
    Δn[end] = (n[end-1] - n[end]) / dz * v_Fend
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
function FE_initalization(bulk_DOS::Union{String, Vector{String}})
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

function athem_magneotransport(sim)
    if sim.athermalelectrons.MagnetoTransport == true
        return :(Δf_mt)
    else
        return :(0.0)
    end
end
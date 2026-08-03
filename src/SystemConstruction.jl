using Base.Threads: @threads
"""
    function_builder(sim::Simulation)
    
    Assembles the correct dictionary of equations for the subsystems that are propagated. Mainly used for
    seperating embedded and un-embedded methods

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Dictionary of subsystem names and their respective equation / expression
"""
function function_builder(sim::Simulation)
    laser=laser_factory(sim)
    if sim.athermalelectrons.EmbeddedAthEM == true
        sim_bg = build_backgroundTTM(sim)
        sys_embed = generate_expressions(sim, laser)
        sys_embed = Dict("$(key)_AthEM" => value for (key, value) in sys_embed)
        sys_bg = generate_expressions(sim_bg, laser)
        comb_sys = merge(sys_embed,sys_bg)
        return comb_sys
    else
        sys = generate_expressions(sim, laser)
        return sys
    end
end
"""
    build_backgroundTTM(sim::Simulation)
    
    Builds the structs for the TTM to be used in the sub-surface layers of the embedding AthEM nmethod

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Simulation struct setup to perform the background TTM for the embedded AthEM method
"""
function build_backgroundTTM(sim::Simulation)
    Tel = ElectronicTemperature(Enabled=true, Electron_PhononCoupling=true, Conductivity=sim.electronictemperature.Conductivity, 
                ElectronicHeatCapacity=sim.electronictemperature.ElectronicHeatCapacity, 
                ElectronPhononCouplingValue = sim.electronictemperature.ElectronPhononCouplingValue)
    Tph = PhononicTemperature(Enabled=true, Electron_PhononCoupling=true, Conductivity=sim.phononictemperature.Conductivity, 
                PhononicHeatCapacity=sim.electronictemperature.PhononicHeatCapacity)
    return build_Simulation(electronictemperature = Tel, phononictemperature=Tph)
end
"""
    generate_expressions(sim::Simulation, laser::Expr)
    
    Calculates and groups each subsystems expression in turn into a dictionary

    # Arguments
    - 'sim': Simulation settings and parameters
    - 'laser': Expression for the temporal evolution and spatial decay of the laser

    # Returns
    - Dictionary of subssytems and their respective expressions
"""
function generate_expressions(sim::Simulation, laser::Expr)
    exprs = Dict{String,Union{Expr,Vector{Expr}}}()
    if sim.electronictemperature.Enabled == true
        merge!(exprs,Dict("Tel" => LightMatter.electrontemperature_factory(sim, laser)))
    end
    if sim.phononictemperature.Enabled == true
        merge!(exprs,Dict("Tph" => LightMatter.phonontemperature_factory(sim)))
    end
    if sim.athermalelectrons.Enabled == true
        merge!(exprs,Dict("fneq" => LightMatter.athemdistribution_factory(sim, laser)))
        if sim.athermalelectrons.AthermalElectron_ElectronCoupling == true
            merge!(exprs,Dict("noe" => LightMatter.athem_thermalelectronparticlechange(sim)))
            τee = electron_relaxationtime(sim)
            merge!(exprs,Dict("relax" => :(LightMatter.athem_electronelectronscattering!(relax_dis, tmp, Tel, μ, sim.structure.egrid, fneq, DOS, noe, $τee, int_vec, Δfexcite))))
        end
    end
    return exprs
end
"""
    generate_initialconditions(sim::Simulation, initialtemps::Dict{String, Float64})
    
    Generates the initial conditions (u0) NamedArrayPartition for the ODE 

    # Arguments
    - 'sim': Simulation settings and parameters
    - 'initialtemps': Dictionary containing initial temepratures for electronic and phononic baths

    # Returns
    - NamedArrayPartition containing the initial conditions of the simulation
"""
function generate_initialconditions(sim::Simulation, initialtemps::Dict{String, Float64})
    return generate_initialconditions(sim, initialtemps, sim.structure)
end

function generate_initialconditions(sim::Simulation, initialtemps::Dict{String, Float64}, ::Structure{1})
    temp_u0 = Dict()
    if sim.athermalelectrons.Enabled == true
        merge!(temp_u0, Dict("fneq" => zeros(sim.structure.dimension.length, length(sim.structure.egrid))))
    end
    if sim.athermalelectrons.AthermalElectron_ElectronCoupling == true
        merge!(temp_u0, Dict("noe" => no_particle_profile(sim)))
    end
    if sim.electronictemperature.Enabled == true
        merge!(temp_u0, Dict("Tel" => fill(initialtemps["Tel"], sim.structure.dimension.length)))
    end
    if sim.phononictemperature.Enabled == true
        merge!(temp_u0, Dict("Tph" => fill(initialtemps["Tph"], sim.structure.dimension.length)))
    end
    namtup = NamedTuple((Symbol(key), value) for (key, value) in temp_u0)
    return NamedArrayPartition(namtup)
end

function generate_initialconditions(sim::Simulation, initialtemps::Dict{String, Float64}, ::Structure{N}) where {N}
    temp_u0 = Dict()
    if sim.athermalelectrons.Enabled == true
        merge!(temp_u0, Dict("fneq" => zeros(sim.structure.dimension.length, length(sim.structure.egrid))))
    end
    if sim.athermalelectrons.AthermalElectron_ElectronCoupling == true
        merge!(temp_u0, Dict("noe" => no_particle_profile(sim)))
    end
    if sim.electronictemperature.Enabled == true
        merge!(temp_u0, Dict("Tel" => fill(initialtemps["Tel"], sim.structure.dimension.length)))
    end
    if sim.phononictemperature.Enabled == true
        merge!(temp_u0, Dict("Tph" => fill(initialtemps["Tph"], sim.structure.dimension.length)))
    end
    namtup = NamedTuple((Symbol(key), value) for (key, value) in temp_u0)
    return NamedArrayPartition(namtup)
end
"""
    generate_initialconditions(sim::Simulation, initialtemps::Dict{String, Float64})
    
    Generates the parameters as a NamedTuple for the ODE 

    # Arguments
    - 'sim': Simulation settings and parameters
    - 'initialtemps': Dictionary containing initial temepratures for electronic and phononic baths

    # Returns
    - NamedTuple containing the parameters of the simulation
"""
function generate_parameters(sim::Simulation, initialtemps::Dict{String, Float64})
    return generate_parameters(sim, initialtemps, sim.structure)
end
###EDITS HERE!
function generate_parameters(sim::Simulation, initialtemps::Dict{String, Float64}, ::Structure{1})
    int_mtx = zeros(sim.structure.dimension.length, length(sim.structure.egrid))
    tmp = zeros(sim.structure.dimension.length, length(sim.structure.egrid))
    p = (sim = sim, int_mtx = int_mtx, tmp = tmp)
    p = parameter_particle(p, sim, sim.structure)
    if sim.athermalelectrons.Enabled
        p = (; p..., Δfexcite = zeros(sim.structure.dimension.length, length(sim.structure.egrid)))
        if !sim.electronictemperature.Enabled
            p = (; p..., Tel = initialtemps["Tel"])
        else
            p = (; p..., relax_dis = zeros(sim.structure.dimension.length, length(sim.structure.egrid)))
        end
    end
    if sim.electronictemperature.Conductivity == true
        p = (; p..., Tel_cond = zeros(sim.structure.dimension.length))
    end
    return p
end

function generate_parameters(sim::Simulation, initialtemps::Dict{String, Float64}, ::Structure{N}) where {N}
    int_mtx = zeros(sim.structure.dimension.length, length(sim.structure.egrid))
    tmp = zeros(sim.structure.dimension.length, length(sim.structure.egrid))
    p = (sim = sim, matsim = sim_seperation(sim), int_mtx = int_mtx, tmp = tmp)
    p = parameter_particle(p, sim, sim.structure)
    if sim.athermalelectrons.Enabled
        p = (; p..., Δfexcite = zeros(sim.structure.dimension.length, length(sim.structure.egrid)))
        if !sim.electronictemperature.Enabled
            p = (; p..., Tel = initialtemps["Tel"])
        else
            p = (; p..., relax_dis = zeros(sim.structure.dimension.length, length(sim.structure.egrid)))
        end
    end
    return p
end

function parameter_particle(p, sim, ::Structure{1})
    if (sim.athermalelectrons.Enabled || sim.electronictemperature.Enabled) && !sim.athermalelectrons.AthermalElectron_ElectronCoupling
        p = (; p..., noe = no_particle_profile(sim))
    end
    return p
end

function parameter_particle(p, sim, ::Structure{N}) where {N}
    if (sim.athermalelectrons.Enabled || sim.electronictemperature.Enabled) && !sim.athermalelectrons.AthermalElectron_ElectronCoupling
        p = (; p..., noe = no_particle_profile(sim))
    end
    return p
end

function no_particle_profile(sim::Simulation)
    if typeof(sim.structure.DOS) == Vector{spl}
        no_part = zeros(sim.structure.dimension.length)
        for j in eachindex(sim.structure.dimension.grid)
            mat = sim.structure.Elemental_System == 1 ? 1 : mat_picker(sim.structure.dimension.grid[j], sim.structure.dimension.InterfaceHeight)
            no_part[j] = get_thermalparticles(0.0, 1e-32, sim.structure.DOS[mat], sim.structure.egrid)
        end
    else
        no_part = fill(get_thermalparticles(0.0, 1e-32, sim.structure.DOS, sim.structure.egrid), sim.structure.dimension.length)
    end
    return no_part
end
"""
    simulation_construction(sys, sim::Simulation)
    
    Creates the expression block for the entire ODE function including multithreading.
    This function calls monometallic_system or antenna_reactor_system respectively.

    # Arguments
    - 'sys': Dictionary of expressions for each subssytem propagated
    - 'sim': Simulation settings and parameters

    # Returns
    - Quote block for the entire ODE problem
"""
function simulation_construction(sys, sim::Simulation, print_time)
    loop_body = build_loopbody(sys, sim)
    if print_time
        t_expr = :(println(t))
    else
        t_expr = :()
    end
    for_header = :(i = 1:p.sim.structure.dimension.length)
    for_expr = Expr(:for, for_header, loop_body)
    if sim.structure.dimension.length > 1
        thread_expr = :(Threads.@threads $for_expr)
    else
        thread_expr = :($for_expr)
    end
    return quote
        $t_expr
        $thread_expr
    end
end
"""
    conductivity_expressions(sim::Simulation)
    
    Creates vector of expressions for the different conductivities occuring during the simulation

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Vector of expression for the conductivity of each subsytem if they are enabled
"""
###EDITS HERE!
function conductivity_expressions(sim::Simulation)
    cond_exprs = [:(LightMatter.reset_du!(du))]
    if sim.electronictemperature.Conductivity == true
        noe_expr = sim.athermalelectrons.Enabled == true && sim.athermalelectrons.AthermalElectron_ElectronCoupling == true ? :(u.noe) : :(p.noe)
        push!(cond_exprs,:(LightMatter.electrontemperature_conductivity!(p.Tel_cond, u.Tel, u.Tph, LightMatter.electronic_heatcapacity_profile(u.Tel, $noe_expr, p.tmp, p.int_mtx, p.sim), p.sim)))
    end
    if sim.phononictemperature.Conductivity == true
        push!(cond_exprs,:(LightMatter.phonontemperature_conductivity!(du.Tph, u.Tph, LightMatter.phononic_heatcapacity_profile(u.Tph, p.sim), p.sim)))
    end
    if sim.athermalelectrons.Conductivity == true && sim.athermalelectrons.AthermalElectron_ElectronCoupling == true
        push!(cond_exprs,:(LightMatter.electron_distribution_transport!(du.fneq, p.sim.athermalelectrons.v_g, u.fneq, p.sim.structure.dimension.spacing, u.Tel, u.noe, p.tmp, p.int_mtx, p.sim)))
    elseif sim.athermalelectrons.Conductivity == true && sim.athermalelectrons.AthermalElectron_ElectronCoupling == false
        push!(cond_exprs,:(LightMatter.electron_distribution_transport!(du.fneq, p.sim.athermalelectrons.v_g, u.fneq, p.sim.structure.dimension.spacing, p.Tel, p.noe, p.tmp, p.int_mtx, p.sim)))
    end
    push!(cond_exprs, :(return nothing))
    return Expr(:block,cond_exprs...)
end

function reset_du!(du)                                                                                                         
    A = ArrayPartition(du)                                                                                                     
    for i in eachindex(A.x)                                                                                                    
        A.x[i] .= 0.0                                                                                                          
    end                                                                                                                        
end
"""
    build_loopbody(sys, sim::Simulation)
    
    Builds the multi-threaded section of the ODE problem

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Quote block for the multithreaded section of the ODE problem
"""
function build_loopbody(sys, sim::Simulation)
    return build_loopbody(sys, sim, sim.structure)
end

function build_loopbody(sys, sim::Simulation, ::Structure{1})
    exprs = Vector{Expr}(undef, 0)
    @debug push!(exprs, :(if ismissing(u) ;  @error "Missing parameter required for simulation." end))
    push!(exprs, variable_renaming(sim))
    if sim.structure.ChemicalPotential
        push!(exprs, :(μ = LightMatter.find_chemicalpotential(noe, Tel, DOS, sim.structure.egrid, tmp, int_vec, μ0)))
    else
        push!(exprs, :(μ = 0.0))
    end
    if sim.athermalelectrons.EmbeddedAthEM == true
        embedding = quote
            if i == 1
                relax_dis = $(sys["relax_AthEM"])
                du.noe = $(sys["noe_AthEM"])
                Δn = du.noe[1]
                du.fneq[1,:] .= $(sys["fneq_AthEM"])
                du.Tel[1] = $(sys["Tel_AthEM"])
                du.Tph[1] = $(sys["Tph_AthEM"])
            else
                du.Tel[i] = $(sys["Tel"])
                du.Tph[i] = $(sys["Tph"])
            end
        end
        push!(exprs,embedding)
    else
        if sim.electronictemperature.Enabled == true && sim.electronictemperature.AthermalElectron_ElectronCoupling == true
            push!(exprs,:($(sys["relax"])))
            push!(exprs,:(@views du.fneq[i,:] .= $(sys["fneq"])))
            push!(exprs,:(@views du.noe[i] = $(sys["noe"])))
        elseif sim.athermalelectrons.Enabled == true
            push!(exprs,:(@views du.fneq[i,:] .= $(sys["fneq"])))
        end
        
        if sim.electronictemperature.Enabled== true
            push!(exprs,:(du.Tel[i] = $(sys["Tel"])))
        end
        if sim.phononictemperature.Enabled == true
            push!(exprs,:(du.Tph[i] = $(sys["Tph"])))
        end
    end
    return Expr(:block, exprs...)
end

function build_loopbody(sys, sim::Simulation, ::Structure{N}) where {N}
    exprs = Vector{Expr}(undef, 0)
    @debug push!(exprs, :(if ismissing(u) ;  @error "Missing parameter required for simulation." end))
    push!(exprs, :(X = LightMatter.mat_picker(p.sim.structure.dimension.grid[i], p.sim.structure.dimension.InterfaceHeight)))
    push!(exprs, ar_variable_renaming(sim))
    if sim.structure.ChemicalPotential
        push!(exprs, :(μ = LightMatter.find_chemicalpotential(noe, Tel, DOS, sim.structure.egrid, tmp, int_vec, μ0)))
    else
        push!(exprs, :(μ = 0.0))
    end
    if sim.athermalelectrons.EmbeddedAthEM == true
        embedding = quote
            if i == 1
                relax_dis = $(sys["relax_AthEM"])
                du.noe = $(sys["noe_AthEM"])
                Δn = du.noe[1]
                du.fneq[1,:] .= $(sys["fneq_AthEM"])
                du.Tel[1] = $(sys["Tel_AthEM"])
                du.Tph[1] = $(sys["Tph_AthEM"])
            else
                du.Tel[i] = $(sys["Tel"])
                du.Tph[i] = $(sys["Tph"])
            end
        end
        push!(exprs, embedding)
    else
        if sim.electronictemperature.Enabled == true && sim.electronictemperature.AthermalElectron_ElectronCoupling == true
            push!(exprs, :($(sys["relax"])))
            push!(exprs, :( @views du.fneq[i,:] .= $(sys["fneq"]) ))
            push!(exprs, :( @views du.noe[i] = $(sys["noe"]) ))
        elseif sim.athermalelectrons.Enabled == true
            push!(exprs, :( @views du.fneq[i,:] .= $(sys["fneq"]) ))
        end

        if sim.electronictemperature.Enabled == true
            push!(exprs, :(du.Tel[i] = $(sys["Tel"])))
        end
        if sim.phononictemperature.Enabled == true
            push!(exprs, :(du.Tph[i] = $(sys["Tph"])))
        end
    end
    return Expr(:block, exprs...)
end
"""
    variable_renaming(sim::Simulation)
    
    Renames variables such as u.x or p.y to just x or y so that during expression evaluation they 
    are named correctly. 

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Expression for the variabble renaming to enter the top of the multithreaded loop
"""
###EDITS HERE!
function variable_renaming(sim::Simulation)
    old_name = [:(p.sim), :(view(p.int_mtx,i,:)), :(@view p.tmp[i,:])]
    new_name = [:sim, :(int_vec), :tmp]
    if typeof(sim.structure.DOS) == Vector{spl}
        push!(old_name, :(p.sim.structure.DOS[i]))
        push!(new_name, :DOS)
    else
        push!(old_name, :(p.sim.structure.DOS))
        push!(new_name, :DOS)
    end
    if sim.athermalelectrons.Enabled == true
        push!(old_name,:(@view u.fneq[i,:]))
        push!(new_name,:fneq)
        push!(old_name, :(@view p.Δfexcite[i,:]))#:(@view LightMatter.access_DiffCache(p.Δfexcite, u.fneq[i,1])[i,:]))
        push!(new_name, :Δfexcite)
        #= if sim.athermalelectrons.Conductivity == true
            push!(old_name, :(@view p.f_cond[i,:]))#:(@view LightMatter.access_DiffCache(p.f_cond, u.fneq[i,1])[i,:]))
            push!(new_name, :f_cond)
        end =#
        if sim.athermalelectrons.AthermalElectron_ElectronCoupling == false
            push!(old_name, :(p.Tel))
            push!(new_name, :Tel)
            push!(old_name, :(p.noe[i]+ LightMatter.get_noparticles(int_vec, fneq, DOS, sim.structure.egrid)))#:(LightMatter.access_DiffCache(p.noe, u.fneq[i,1])[i]))
            push!(new_name, :n)
        else 
            push!(old_name, :(u.noe[i] + LightMatter.get_noparticles(int_vec, fneq, DOS, sim.structure.egrid)))
            push!(new_name, :noe)
            push!(old_name, :(@view p.relax_dis[i,:]))#:(@view LightMatter.access_DiffCache(p.relax_dis, u.fneq[i,1])[i,:]))
            push!(new_name, :relax_dis)
        end
    end
    if sim.electronictemperature.Enabled == true
        push!(old_name, :(u.Tel[i]))
        push!(new_name, :Tel)
        if sim.athermalelectrons.AthermalElectron_ElectronCoupling == false
            push!(old_name, :(p.noe[i]))#:(LightMatter.access_DiffCache(p.noe, u.Tel[i])[i]))
            push!(new_name, :noe)
        end
        if sim.electronictemperature.Conductivity == true
            push!(old_name,:(p.Tel_cond[i]))#:(LightMatter.access_DiffCache(p.Tel_cond,u.Tel[i])[i]))
            push!(new_name,:Tel_cond)
        end
    end
    if sim.phononictemperature.Enabled == true
        push!(old_name, :(u.Tph[i]))
        push!(new_name, :Tph)
        #= if sim.phononictemperature.Conductivity == true
            push!(old_name, :(p.Tph_cond[i]))#:(LightMatter.access_DiffCache(p.Tph_cond,u.Tph[i])[i]))
            push!(new_name, :Tph_cond)
        end =#
    end
    old_name = Tuple(old_name)
    new_name = Tuple(new_name)
    assignments = [:( local $(lhs) = $(rhs)) for (lhs, rhs) in zip(new_name, old_name)]
    return quote
        $(assignments...)
    end
end
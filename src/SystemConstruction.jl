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
            merge!(exprs,Dict("relax" => :(LightMatter.athem_electronelectronscattering!(relax_dis, tmp, Tel, μ, sim.structure.egrid, fneq, DOS, n, $τee))))
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
    temp_u0 = Dict()
    if sim.athermalelectrons.Enabled == true
        merge!(temp_u0, Dict("fneq"=>zeros(sim.structure.dimension.length, length(sim.structure.egrid))))
    end
    if sim.athermalelectrons.AthermalElectron_ElectronCoupling == true 
        if typeof(sim.structure.DOS) == Vector{spl}
            if sim.structure.Elemental_System == 1
                no_part=zeros(sim.structure.dimension.length)
                for j in eachindex(sim.structure.dimension.grid)
                    no_part[j] = get_thermalparticles(0.0, 1e-32, sim.structure.DOS[j], sim.structure.egrid)
                end
            else
                no_part=zeros(sim.structure.dimension.length)
                for j in eachindex(sim.structure.dimension.grid)
                    x = mat_picker(sim.structure.dimension.grid[j], sim.structure.dimension.InterfaceHeight)
                    no_part[j] = get_thermalparticles(0.0, 1e-32, sim.structure.DOS[x], sim.structure.egrid)
                end
            end
        else
            no_part = fill(get_thermalparticles(0.0,1e-32, sim.structure.DOS, sim.structure.egrid), sim.structure.dimension.length)
        end
        merge!(temp_u0,Dict("noe" => no_part ))
    end
    if sim.electronictemperature.Enabled == true
        merge!(temp_u0, Dict("Tel" => fill(initialtemps["Tel"], sim.structure.dimension.length)))
    end
    if sim.phononictemperature.Enabled == true
        merge!(temp_u0, Dict("Tph" => fill(initialtemps["Tph"], sim.structure.dimension.length)))
    end
    namtup = NamedTuple((Symbol(key),value) for (key,value) in temp_u0)
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
    if sim.structure.Elemental_System > 1
        p = (sim=sim,matsim=sim_seperation(sim))
    else
        p = (sim=sim,)
    end
    p = parameter_conductivity(p, sim)
    p = parameter_particle(p, sim)
    if sim.athermalelectrons.Enabled
        p =(; p..., Δfexcite = [DiffCache(zeros(length(sim.structure.egrid))) for _ in 1:sim.structure.dimension.length],
                    tmp = [DiffCache(zeros(length(sim.structure.egrid))) for _ in 1:sim.structure.dimension.length])
        if !sim.electronictemperature.Enabled
            p = (; p..., Tel = initialtemps["Tel"])
        else
            p = (; p..., relax_dis = [DiffCache(zeros(length(sim.structure.egrid))) for _ in 1:sim.structure.dimension.length])
        end
    end

    return p 
end

function parameter_conductivity(p, sim)
    if sim.electronictemperature.Conductivity
        p = (; p..., Tel_cond = DiffCache(zeros(sim.structure.dimension.length)))
    end
    if sim.athermalelectrons.Conductivity
        p = (; p..., f_cond = DiffCache(zeros(sim.structure.dimension.length, length(sim.structure.egrid))))
    end
    if sim.phononictemperature.Conductivity
        p = (; p..., Tph_cond = DiffCache(zeros(sim.structure.dimension.length)))
    end
    return p
end

function parameter_particle(p, sim)
    if (sim.athermalelectrons.Enabled || sim.electronictemperature.Enabled) && !sim.athermalelectrons.AthermalElectron_ElectronCoupling
        if typeof(sim.structure.DOS) == Vector{spl} && sim.structure.Elemental_System == 1
            no_part=zeros(sim.structure.dimension.length)
            noe = get_tmp(no_part, 0.0)
            for j in eachindex(sim.structure.dimension.grid)
                noe[j] = get_thermalparticles(0.0, 1e-32, sim.structure.DOS[j], sim.structure.egrid)
            end
        elseif typeof(sim.structure.DOS) == Vector{spl} && sim.structure.Elemental_System > 1
            no_part=DiffCache(zeros(sim.structure.dimension.length))
            noe = get_tmp(no_part, 0.0)
            for j in eachindex(sim.structure.dimension.grid)
                mat = mat_picker(sim.structure.dimension.grid[j], sim.structure.dimension.InterfaceHeight)
                noe[j] = get_thermalparticles(0.0, 1e-32, sim.structure.DOS[mat], sim.structure.egrid)
            end
        else
            no_part = DiffCache(fill(get_thermalparticles(0.0, 1e-32, sim.structure.DOS, sim.structure.egrid), sim.structure.dimension.length))
        end
        p = (; p..., noe = no_part)
    end
    return p
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
    if sim.structure.Elemental_System == 1
        return monometallic_system(sys, sim, print_time)
    else
        return antenna_reactor_system(sys, sim, print_time)
    end
end
"""
    monometallic_system(sys, sim::Simulation)
    
    Creates the expression block for the entire ODE function including multithreading.
    This function is sepcfically for a monometallic system, Elemental_System == 1

    # Arguments
    - 'sys': Dictionary of expressions for each subssytem propagated
    - 'sim': Simulation settings and parameters

    # Returns
    - Quote block for the monometallic problem
"""
function monometallic_system(sys, sim::Simulation, print_time)
    expr_cond = conductivity_expressions(sim)
    loop_body = build_loopbody(sys, sim)
    if print_time
        t_expr = :(println(t))
    else
        t_expr = :()
    end
    return quote
        $t_expr
        $expr_cond
        Threads.@threads for i in 1:p.sim.structure.dimension.length
            $loop_body
        end
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
function conductivity_expressions(sim::Simulation)
    cond_exprs = []
    if sim.electronictemperature.Conductivity == true
        push!(cond_exprs,:(LightMatter.electrontemperature_conductivity!(u.Tel, p.sim.electronictemperature.κ, p.sim.structure.dimension.spacing, u.Tph, p.Tel_cond)))
    end
    if sim.phononictemperature.Conductivity == true
        push!(cond_exprs,:(LightMatter.phonontemperature_conductivity!(u.Tph, p.sim.phononictemperature.κ, p.sim.structure.dimension.spacing, p.Tph_cond)))
    end
    if sim.athermalelectrons.Conductivity == true
        push!(cond_exprs,:(LightMatter.electron_distribution_transport!(p.sim.athermalelectrons.v_g, u.fneq, p.f_cond, p.sim.structure.dimension.spacing)))
    end 
    return Expr(:block,cond_exprs...)
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
    exprs = Vector{Expr}(undef,0)
    if sim.structure.Elemental_System > 1
        push!(exprs, :(X = LightMatter.mat_picker(p.sim.structure.dimension.grid[i], p.sim.structure.dimension.InterfaceHeight))) # Picks the active material
        push!(exprs, ar_variable_renaming(sim)) # Translates variable names from DiffEq.jl to LightMatter.jl
    else
        push!(exprs,variable_renaming(sim))
    end
    if sim.structure.ChemicalPotential
        push!(exprs, :(μ = LightMatter.find_chemicalpotential(n, Tel, DOS, sim.structure.egrid)))
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
"""
    variable_renaming(sim::Simulation)
    
    Renames variables such as u.x or p.y to just x or y so that during expression evaluation they 
    are named correctly. 

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Expression for the variabble renaming to enter the top of the multithreaded loop
"""
function variable_renaming(sim::Simulation)
    old_name = [:(p.sim)]
    new_name = [:sim]
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
        push!(old_name, :(p.tmp[i]))
        push!(new_name, :tmp)
        push!(old_name, :(p.Δfexcite[i]))
        push!(new_name, :Δfexcite)
        if sim.athermalelectrons.Conductivity == true
            push!(old_name,:(@view LightMatter.access_DiffCache(p.f_cond, u.fneq[i,1])[i,:]))
            push!(new_name,:f_cond)
        end
        if sim.athermalelectrons.AthermalElectron_ElectronCoupling == false
            push!(old_name,:(p.Tel))
            push!(new_name,:Tel)
            push!(old_name, :(LightMatter.access_DiffCache(p.noe, u.fneq[i,1])[i]))
            push!(new_name, :n)
        else 
            push!(old_name, :(u.noe[i] + LightMatter.get_noparticles(fneq, DOS, sim.structure.egrid)))
            push!(new_name, :n)
            push!(old_name, :(p.relax_dis[i]))
            push!(new_name, :relax_dis)
        end
    end
    if sim.electronictemperature.Enabled == true
        push!(old_name,:(u.Tel[i]))
        push!(new_name,:Tel)
        if sim.athermalelectrons.AthermalElectron_ElectronCoupling == false
            push!(old_name, :(LightMatter.access_DiffCache(p.noe, u.Tel[i])[i]))
            push!(new_name, :n)
        end
        if sim.electronictemperature.Conductivity == true
            push!(old_name,:(LightMatter.access_DiffCache(p.Tel_cond,u.Tel[i])[i]))
            push!(new_name,:Tel_cond)
        end
    end
    if sim.phononictemperature.Enabled == true
        push!(old_name,:(u.Tph[i]))
        push!(new_name,:Tph)
        if sim.phononictemperature.Conductivity == true
            push!(old_name,:(LightMatter.access_DiffCache(p.Tph_cond,u.Tph[i])[i]))
            push!(new_name,:Tph_cond)
        end
    end
    old_name = Tuple(old_name)
    new_name = Tuple(new_name)
    assignments = [:( local $(lhs) = $(rhs)) for (lhs, rhs) in zip(new_name, old_name)]
    return quote
        $(assignments...)
    end
end
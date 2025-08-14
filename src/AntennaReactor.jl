"""
    antenna_reactor_system(sys, sim)

    Constructs the main expression for an antenna reactor simulation system.

    # Arguments
    - `sys`::Dict{String, Union{Expr, Vector{Expr}}}
           A dictionary or configuration object containing the expressions for each system propagated
    - `sim`::Simulation
           A simulation object containing all necessary input data for the simulation.

    # Returns
    - An expression block containing the time output, conductivity expressions, and a threaded simulation loop.
"""
function antenna_reactor_system(sys::Dict{String, Union{Expr, Vector{Expr}}}, sim::Simulation, print_time)
    expr_cond = conductivity_expressions(sim) # Expr block for the conductivity of each system
    loop_body = build_loopbody(sys, sim) # Expr block for the body of a threaded for loop over the systems and depth
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
    end # The whole simulation that is propgated as one Expr block
end
"""
    mat_picker(height, cutoffs)

    Selects an index based on material interface height and given cutoffs.

    # Arguments
    - `height`::Float64
              The height value to compare.
    - `cutoffs`::Union{Float64,Vector{Float64}}
               A vector of cutoff heights defining material regions.

    # Returns
    - The index of the region in which `height` lies.
"""
function mat_picker(height::Float64, cutoffs::Union{Float64,Vector{Float64}})
    subindex = 1
    for i in eachindex(cutoffs) 
        if height < cutoffs[i]
            subindex = i
            break
        else
            subindex = length(cutoffs) + 1  # If beyond the highest cutoff
        end
    end
    return subindex
end
"""
    ar_variable_renaming(sim)

    Generates variable renaming expressions for translation between the variibale names in DiffEq.jl and LightMatter.jl

    # Arguments
    - `sim`::Simulation
           The Simulation struct containing all information about the simulation

    # Returns
    - An expression block assigning simulation-specific variable names.
"""
function ar_variable_renaming(sim::Simulation)
    old_name = [:(p.matsim[X])]
    new_name = [:sim]
    if typeof(sim.structure.DOS) == Vector{Vector{spl}}
        push!(old_name, :(p.matsim[X].structure.DOS[i]))
        push!(new_name, :DOS)
    else
        push!(old_name, :(p.matsim[X].structure.DOS))
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
            push!(old_name, :(@view LightMatter.access_DiffCache(p.f_cond, u.fneq[i,1])[i,:]))
            push!(new_name, :f_cond)
        end
        if sim.athermalelectrons.AthermalElectron_ElectronCoupling == false
            push!(old_name, :(p.Tel))
            push!(new_name, :Tel)
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
        push!(old_name, :(u.Tel[i]))
        push!(new_name, :Tel)
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
        push!(old_name, :(u.Tph[i]))
        push!(new_name, :Tph)
        if sim.phononictemperature.Conductivity == true
            push!(old_name, :(LightMatter.access_DiffCache(p.Tph_cond,u.Tph[i])[i]))
            push!(new_name, :Tph_cond)
        end
    end
    old_name = Tuple(old_name)
    new_name = Tuple(new_name)
    assignments = [:(local $(lhs) = $(rhs)) for (lhs, rhs) in zip(new_name, old_name)]
    return quote
        $(assignments...)
    end
end
"""
    sim_seperation(sim)

    Splits a composite `Simulation` object into separate simulations for each elemental subsystem.

    # Arguments
    - `sim`::Simulation
           The composite simulation containing multiple subsystems.

    # Returns
    - A vector of `Simulation` objects, each corresponding to a single elemental subsystem.
"""
function sim_seperation(sim::Simulation)
    function conditional_split(field, enabled::Bool) #Splits the system into its subcomponents or replicates the disabled system
        enabled ? split_struct(field, sim.structure.Elemental_System) :
                  fill(field, sim.structure.Elemental_System)
    end
    lasers = fill(sim.laser, sim.structure.Elemental_System) #Will always be enabled so doesn't need checking
    tels   = conditional_split(sim.electronictemperature, sim.electronictemperature.Enabled)
    tphs   = conditional_split(sim.phononictemperature, sim.phononictemperature.Enabled)
    neqs   = conditional_split(sim.athermalelectrons, sim.athermalelectrons.Enabled)
    structs = split_structure(sim.structure)
    new_sim = Vector{SimulationTypes}(undef, sim.structure.Elemental_System) #Build and fill vector of simulation objects
    for i in 1:sim.structure.Elemental_System 
        new_sim[i] = build_Simulation(laser=lasers[i], electronictemperature=tels[i], phononictemperature=tphs[i],
                                    athermalelectrons=neqs[i], structure=structs[i])
    end
    return new_sim
end
"""
    split_struct(data, number)

    Splits fields of a composite object with vector fields into a vector of scalar instances.

    # Arguments
    - `data`::SimulationTypes 
            A subsystem of the simulation e.g. AthermalElectrons or ElectronicTemperature
    - `number`::Int 
              Number of different elemental systems.

    # Returns
    - A vector of the subsystem with scalar data extracted from the original vector fields.
"""
function split_struct(data::SimulationTypes, number::Int)
    field_values = map(f -> getfield(data, f), fieldnames(typeof(data))) #Extracts all values from fields of subssytem into dict
    value_indices = findall(x -> x isa Vector, field_values) # Determines which fields of the subsystem need splitting 

    return [
        typeof(data)(
            map((val, idx) -> idx in value_indices && length(val) != 1 ? val[i] : val, field_values, 1:length(field_values))...
        ) for i in 1:number
    ] # Creates a vector of the subsystem with each containing one value of the vector fields as a scalar
end
"""
    split_structure(structure)

    Splits a `Structure` object into multiple structures if `DOS` is a vector. Different to split_struct due to the possibility 
    of spatially resolved DOS' though currently that isn't implemented

    # Arguments
    - `structure`::Structure
                 The struct containing the information of the simulatkion structure

    # Returns
    - A vector of `Structure` objects, one per element if applicable.
"""
function split_structure(structure::Structure)
    field_values = Dict(f => getfield(structure, f) for f in fieldnames(typeof(structure))) #Extracts all values from fields of subssytem into dict
    
    if !(:DOS in keys(field_values) && field_values[:DOS] isa Vector)
        return fill(structure, structure.Elemental_System)  # Return as-is if the specified field is not a vector
    end
    
    return [
        typeof(structure)(
            (f == :DOS ? field_values[f][i] : field_values[f] for f in fieldnames(typeof(structure)))...
        ) for i in 1:length(field_values[:DOS])
    ] # Creates a vector of the Structure struct with the DOS split into their seperate materials
end
"""
    split_grid(grid, cutoffs)

    Splits a numerical grid into regions that connect to each material

    # Arguments
    - `grid`::Vector{Float64}
            A vector of Float64 numbers representing the full z-grid
    - `cutoffs`::Union{Float64, Vector{Float64}}
               A single value or vector of values defining the interfaces between each material 

    # Returns
    - A vector of sub-vectors representing segments of the original z-grid.
"""
function split_grid(grid::Vector{Float64}, cutoffs::Union{Float64,Vector{Float64}})
    cutoffs = isa(cutoffs, Vector) ? cutoffs : [cutoffs]
    sections = Vector{Vector{Float64}}()
    start_idx = 1
    for cutoff in cutoffs
        end_idx = findfirst(x -> x > cutoff, grid)
        if isnothing(end_idx)
            push!(sections, grid[start_idx:end])
            return sections
        end
        push!(sections, grid[start_idx:end_idx-1])
        start_idx = end_idx
    end
    push!(sections, grid[start_idx:end])
    return sections
end

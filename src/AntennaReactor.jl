###
    #To-Do: Add spatial DOS support to antenna-reactor complex
###
"""
    antenna_reactor_system(sys::Dict{String, Union{Expr, Vector{Expr}}}, sim::Simulation)

    Constructs the main expression for an antenna reactor simulation system.

    # Arguments
    - `sys`: A dictionary or configuration object containing the expressions for each system propagated
    - `sim`: A simulation object containing all necessary input data for the simulation.

    # Returns
    - An expression block containing the time output, conductivity expressions, and a threaded simulation loop.
"""
function antenna_reactor_system(sys::Dict{String, Union{Expr, Vector{Expr}}}, sim::Simulation)
    expr_cond = conductivity_expressions(sim) # Expr block for the conductivity of each system
    loop_body = ar_build_loopbody(sys, sim) # Expr block for the body of a threaded for loop over the systems and depth
    return quote 
        println(t)
        $expr_cond
        Threads.@threads for i in 1:p.sim.structure.dimension.length
            $loop_body
        end
    end # The whole simulation that is propgated as one Expr block
end
"""
    ar_build_loopbody(sys::Dict{String, Union{Expr, Vector{Expr}}}, sim::Simulation)

    Builds the core loop body expression for the antenna reactor simulation.

    # Arguments
    - `sys`: Dict of Expr for each subsystem
    - `sim`: All simulation parameters

    # Returns
    - An expression representing the simulation loop body.
"""
function ar_build_loopbody(sys::Dict{String, Union{Expr, Vector{Expr}}}, sim::Simulation)
    exprs = Vector{Expr}(undef, 0) # Miscellaneous expressions at the top of the loop
    push!(exprs, :(X = Lightmatter.mat_picker(p.sim.structure.dimension.grid[i], p.sim.structure.dimension.InterfaceHeight))) # Picks the active material
    push!(exprs, ar_variable_renaming(sim)) # Translates variable names from DiffEq.jl to Lightmatter.jl
    push!(exprs, :(μ = Lightmatter.find_chemicalpotential(n, Tel, DOS, sim.structure.egrid))) # Calculates the current chemical potential

    if sim.athermalelectrons.EmbeddedAthEM == true 
        # Embeds AthEM to reduce computational cost - do not use with non-eq transport but with both electronic and phononic temperature
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
    else # No embedding so all heights are treated the same 
        if sim.electronictemperature.Enabled == true && sim.electronictemperature.AthermalElectron_ElectronCoupling == true
            push!(exprs,:(relax_dis = $(sys["relax"])))
            push!(exprs,:(du.noe[i,:] .= $(sys["noe"])))
            push!(exprs,:(Δn = du.noe[i]))
        end

        if sim.athermalelectrons.Enabled == true
            push!(exprs,:(du.fneq[i,:] .= $(sys["fneq"])))
        end
        if sim.electronictemperature.Enabled== true
            push!(exprs,:(du.Tel[i,:] .= $(sys["Tel"])))
        end
        if sim.phononictemperature.Enabled == true
            push!(exprs,:(du.Tph[i,:] .= $(sys["Tph"])))
        end
    end
    return Expr(:block, exprs...)
end
"""
    mat_picker(height::Real, cutoffs::Union{Real,Vector{<:Real}})

    Selects an index based on material interface height and given cutoffs.

    # Arguments
    - `height`: The height value to compare.
    - `cutoffs`: A vector of cutoff heights defining material regions.

    # Returns
    - The index of the region in which `height` lies.
"""
function mat_picker(height::Real, cutoffs::Union{Real,Vector{<:Real}})
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
    ar_variable_renaming(sim::Simulation)

    Generates variable renaming expressions for translation between the variibale names in DiffEq.jl and Lightmatter.jl

    # Arguments
    - `sim`: The Simulation struct containing all information about the simulation

    # Returns
    - An expression block assigning simulation-specific variable names.
"""
function ar_variable_renaming(sim::Simulation)
    old_name = [:(p.matsim[X])]
    new_name = [:sim,]
    if typeof(sim.structure.DOS) == Vector{Vector{spl}}
        push!(old_name, :(p.matsim[X].structure.DOS[i]))
        push!(new_name, :DOS)
    else
        push!(old_name, :(p.matsim[X].structure.DOS))
        push!(new_name, :DOS)
    end
    if sim.athermalelectrons.Enabled == true
        push!(old_name, :(u.fneq[i,:]))
        push!(new_name, :fneq)
        if sim.athermalelectrons.Conductivity == true
            push!(old_name, :(p.f_cond[i,:]))
            push!(new_name, :f_cond)
        end
        if sim.athermalelectrons.AthermalElectron_ElectronCoupling == false
            push!(old_name, :(p.Tel))
            push!(new_name, :Tel)
            push!(old_name, :(p.noe[i]))
            push!(new_name, :n)
        else 
            push!(old_name, :(u.noe[i]))
            push!(new_name, :n)
        end
    end
    if sim.electronictemperature.Enabled == true
        push!(old_name, :(u.Tel[i]))
        push!(new_name, :Tel)
        if sim.athermalelectrons.AthermalElectron_ElectronCoupling == false
            push!(old_name, :(p.noe[i]))
            push!(new_name, :n)
        end
        if sim.electronictemperature.Conductivity == true
            push!(old_name,:(p.Tel_cond[i]))
            push!(new_name,:Tel_cond)
        end
    end
    if sim.phononictemperature.Enabled == true
        push!(old_name, :(u.Tph[i]))
        push!(new_name, :Tph)
        if sim.phononictemperature.Conductivity == true
            push!(old_name, :(p.Tph_cond[i]))
            push!(new_name, :Tph_cond)
        end
    end
    old_name = Tuple(old_name)
    new_name = Tuple(new_name)
    assignments = [:(local $(lhs) = $(rhs)) for (lhs, rhs) in zip(new_name, old_name)]
    return quote
        $(assignments...)
    end # Writes a single line which renames all the variables correctly
end
"""
    sim_seperation(sim::Simulation)

    Splits a composite `Simulation` object into separate simulations for each elemental subsystem.

    # Arguments
    - `sim`: The composite simulation containing multiple subsystems.

    # Returns
    - A vector of `Simulation` objects, each corresponding to a single elemental subsystem.
"""
function sim_seperation(sim::Simulation)
    function conditional_split(field, enabled::Bool) #Splits the system into its subcomponents or replicates the disabled system
        enabled ? split_struct(field, sim.structure.Elemental_System) :
                  fill(field, sim.structure.Elemental_System)
    end
    lasers = split_struct(sim.laser, sim.structure.Elemental_System) #Will always be enabled so doesn't need checking
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
    split_struct(data::SimulationTypes, number::Int)

    Splits fields of a composite object with vector fields into a vector of scalar instances.

    # Arguments
    - `data`: A subsystem of the simulation e.g. neq electrons or Tel
    - `number`: Number of subdivisions (usually the number of elements).

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
    split_structure(structure::Structure)

    Splits a `Structure` object into multiple structures if `DOS` is a vector. Different to split_struct due to the possibility 
    of spatially resolved DOS' though currently that isn't implemented

    # Arguments
    - `structure`: The struct containing the information of the simulatkion structure

    # Returns
    - A vector of `Structure` objects, one per element if applicable.
"""
function split_structure(structure::Structure)
    field_values = Dict(f => getfield(structure, f) for f in fieldnames(typeof(structure))) #Extracts all values from fields of subssytem into dict
    
    if !(:DOS in keys(field_values) && field_values[:DOS] isa Vector)
        return [structure]  # Return as-is if the specified field is not a vector
    end
    
    return [
        typeof(structure)(
            (f == :DOS ? field_values[f][i] : field_values[f] for f in fieldnames(typeof(structure)))...
        ) for i in 1:length(field_values[:DOS])
    ] # Creates a vector of the Structure struct with the DOS split into their seperate materials
end
"""
    split_grid(grid::Vector{<:Real}, cutoffs::Union{Real, Vector{Real}})

    Splits a numerical grid into regions that connect to each material

    # Arguments
    - `grid`: A vector of real numbers representing the full z-grid
    - `cutoffs`: A single value or vector of values defining the interfaces between each material 

    # Returns
    - A vector of sub-vectors representing segments of the original z-grid.
"""
function split_grid(grid::Vector{<:Real}, cutoffs::Union{Real,Vector{Real}})
    cutoffs = isa(cutoffs, Vector) ? cutoffs : [cutoffs]
    sections = Vector{Vector{<:Real}}()
    start_idx = 1
    for cutoff in sorted_cutoffs
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

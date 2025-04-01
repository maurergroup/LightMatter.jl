function antenna_reactor_system(sys,sim::Simulation)
    cond_exprs = conductivity_expressions(sim)
    loop_body = build_loopbody(sys,sim)
    expr_cond = Expr(:block,cond_exprs...)
    return quote
        println(t)
        $expr_cond
        Threads.@threads for i in 1:p.sim.structure.dimension.length
            $loop_body
        end
    end
end

function sim_seperation(sim::Simulation)
    lasers = split_struct(sim.laser)
    tels = split_struct(sim.electronictemperature)
    tphs = split_struct(sim.phononictemeprature)
    neqs = split_struct(sim.athermalelectrons)
    structs = split_structure(sim.structure)

    new_sim = Vector{Simulation}(undef,1:sim.structure.Elemental_System)
    for i in 1:sim.structure.Elemental_System
        new_sim[i] = build_Simulation(laser=lasers[i],electronictemperature=tels[i],phononictemperature=tphs[i],
                                    athermalelectrons=neqs[i],structure=structs[i])
    end
    return new_sim
end

function split_struct(data::SimulationTypes)
    field_values = map(f -> getfield(data, f), fieldnames(typeof(data)))
    value_indices = findall(x -> x isa Vector, field_values)
    
    if isempty(value_indices)
        return [data]  # Return as-is if no vector fields exist
    end
    
    vec_lengths = map(i -> length(field_values[i]), value_indices)
    n = maximum(vec_lengths)
    
    return [
        typeof(data)(
            map((val, idx) -> idx in value_indices ? [val[i]] : val, field_values, 1:length(field_values))...
        ) for i in 1:n
    ]
end

function split_structure(structure::Structure)
    field_values = Dict(f => getfield(structure, f) for f in fieldnames(typeof(structure)))
    
    if !(:DOS in keys(field_values) && field_values[:DOS] isa Vector)
        return [structure]  # Return as-is if the specified field is not a vector
    end
    
    return [
        typeof(structure)(
            (f == :DOS ? [field_values[f][i]] : field_values[f] for f in fieldnames(typeof(structure)))...
        ) for i in 1:length(field_values[:DOS])
    ]
end

function cut_offloop()
    for h in heights
        # Determine which region the current height falls into
        for i in 1:length(cutoffs)
            if h < cutoffs[i]
                subindex = i
                break
            else
                subindex = length(cutoffs) + 1  # If beyond the highest cutoff
            end
        end
        println("Height: $h, Subindex: $subindex")
    end
end
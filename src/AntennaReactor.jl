function antenna_reactor_system(sys,sim::Simulation)
    cond_exprs = conductivity_expressions(sim)
    loop_body = ar_build_loopbody(sys,sim)
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
    lasers = split_struct(sim.laser,sim.structure.Elemental_System)
    if sim.electronictemperature.Enabled == true
        tels = split_struct(sim.electronictemperature,sim.structure.Elemental_System)
    else
        tels = fill(sim.electronictemperature,sim.structure.Elemental_System)
    end
    if sim.phononictemperature.Enabled == true
        tphs = split_struct(sim.phononictemperature,sim.structure.Elemental_System)
    else
        tphs = fill(sim.phononictemperature,sim.structure.Elemental_System)
    end
    if sim.athermalelectrons.Enabled == true
        neqs = split_struct(sim.athermalelectrons,sim.structure.Elemental_System)
    else
        neqs = fill(sim.athermalelectrons,sim.structure.Elemental_System)
    end
    structs = split_structure(sim.structure)

    new_sim = Vector{SimulationTypes}(undef,sim.structure.Elemental_System)
    for i in 1:sim.structure.Elemental_System
        new_sim[i] = build_Simulation(laser=lasers[i],electronictemperature=tels[i],phononictemperature=tphs[i],
                                    athermalelectrons=neqs[i],structure=structs[i])
    end
    return new_sim
end

function split_struct(data::SimulationTypes,number)
    field_values = map(f -> getfield(data, f), fieldnames(typeof(data)))
    value_indices = findall(x -> x isa Vector, field_values)
    
    if isempty(value_indices)
        return [data]  # Return as-is if no vector fields exist
    end

    return [
        typeof(data)(
            map((val, idx) -> idx in value_indices && length(val) != 1 ? val[i] : val, field_values, 1:length(field_values))...
        ) for i in 1:number
    ]
end

function split_structure(structure::Structure)
    field_values = Dict(f => getfield(structure, f) for f in fieldnames(typeof(structure)))
    
    if !(:DOS in keys(field_values) && field_values[:DOS] isa Vector)
        return [structure]  # Return as-is if the specified field is not a vector
    end
    
    return [
        typeof(structure)(
            (f == :DOS ? field_values[f][i] : field_values[f] for f in fieldnames(typeof(structure)))...
        ) for i in 1:length(field_values[:DOS])
    ]
end

function mat_picker(height, cutoffs)
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

function ar_build_loopbody(sys,sim::Simulation)
    exprs = Vector{Expr}(undef,0)
    push!(exprs,:(X = Lightmatter.mat_picker(p.sim.structure.dimension.grid[i],p.sim.structure.dimension.InterfaceHeight)))
    push!(exprs,ar_variable_renaming(sim))
    push!(exprs, :(μ = Lightmatter.find_chemicalpotential(n,Tel,DOS,sim.structure.egrid)))

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
        push!(old_name,:(u.fneq[i,:]))
        push!(new_name,:fneq)
        if sim.athermalelectrons.Conductivity == true
            push!(old_name,:(p.f_cond[i,:]))
            push!(new_name,:f_cond)
        end
        if sim.athermalelectrons.AthermalElectron_ElectronCoupling == false
            push!(old_name,:(p.Tel))
            push!(new_name,:Tel)
            push!(old_name, :(p.noe[i]))
            push!(new_name, :n)
        else 
            push!(old_name, :(u.noe[i]))
            push!(new_name, :n)
        end
    end
    if sim.electronictemperature.Enabled == true
        push!(old_name,:(u.Tel[i]))
        push!(new_name,:Tel)
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
        push!(old_name,:(u.Tph[i]))
        push!(new_name,:Tph)
        if sim.phononictemperature.Conductivity == true
            push!(old_name,:(p.Tph_cond[i]))
            push!(new_name,:Tph_cond)
        end
    end
    old_name = Tuple(old_name)
    new_name = Tuple(new_name)
    assignments = [:(local $(lhs) = $(rhs)) for (lhs, rhs) in zip(new_name, old_name)]
    return quote
        $(assignments...)
    end
end
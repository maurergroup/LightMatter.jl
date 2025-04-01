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
end

function split_struct(data)
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

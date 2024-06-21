using Symbolics,ModelingToolkit

function summing_vectors(vec1::Vector{Float64},vec2::Vector{Float64})
    return vec1.+vec2
end
@register_array_symbolic summing_vectors(vec1::AbstractVector,vec2::AbstractVector) begin
    size = (length(vec1),)
    eltype= eltype(vec1)
end

@variables (output(t))[1:200] (input_1(t))[1:200] (input_2(t))[1:200]

output = summing_vectors(input_1,input_2)
size(output)

#= subs = Dict([input_1 => range(0,10,length=200),input_2 => range(20,40,length=200)])

answer = substitute(eq,subs) =#
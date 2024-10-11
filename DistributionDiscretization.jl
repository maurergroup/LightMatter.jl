using Random,Distributions

function generate_DOS(File::String,n)
    TotalDOS::Matrix{Float64}=readdlm(File)
    return get_interpolate(TotalDOS[:,1],TotalDOS[:,2].*n)
end

function split_into_binary_vectors(x::Vector{Float64}, n::Int, w::Vector{Float64})
    # Ensure x is between 0 and 1
    @assert all(0 .<= x .<= 1) "All elements of x must be between 0 and 1"
    @assert length(x) == length(w) "x and w must have the same length"
    
    m = length(x)  # Length of the input vector
    binary_vectors = zeros(Int, m, n)  # Initialize binary vectors

    # Step 1: Generate binary vectors using probabilistic rounding
    for i in 1:m
        for j in 1:n
            binary_vectors[i, j] = rand(Bernoulli(x[i]))  # Random binary based on probability x[i]
        end
    end

    # Step 2: Compute the original weighted sum
    original_weighted_sum = sum(x .* w)

    # Step 3: Adjust each binary vector to match the weighted sums of the original vector
    for j in n:1
        # Compute the current weighted sum for the j-th binary vector
        current_weighted_sum = sum(binary_vectors[:, j] .* w)
        diff = original_weighted_sum - current_weighted_sum

        # Step 4: Adjust the binary vector to correct for the difference
        for i in 1:m
            if diff == 0
                break
            elseif diff > 0 && binary_vectors[i, j] == 0 && w[i] > 0
                # Flip 0 -> 1 if it helps reduce the positive difference
                binary_vectors[i, j] = 1
                diff -= w[i]
            elseif diff < 0 && binary_vectors[i, j] == 1 && w[i] > 0
                # Flip 1 -> 0 if it helps reduce the negative difference
                binary_vectors[i, j] = 0
                diff += w[i]
            end
        end
    end

    # Step 5: Return the resulting binary vectors
    return binary_vectors
end

function restore_ftot(splits,n)
    restore=zeros(length(splits[:,1]))
    for i in eachindex(splits[:,1])
        restore[i]=sum(splits[i,:])
    end
    return restore./n
end

egrid = collect(range(-3*1.55,3*1.55,step=0.02))
DOS_spl = generate_DOS("DOS/Cu_DOS.dat",59)
DOS = DOS_spl(egrid)
fneq = sol[4:end,451]
Tel=sol[1,451]
n=sol[2,451]
μ=find_chemicalpotential(n,Tel,DOS_spl,cons.kB,mp.FE,mp.n0)
feq = FermiDirac(Tel,μ,cons.kB,egrid)
ftot = feq.+fneq

vecs=10000
splits = split_into_binary_vectors(ftot,vecs,DOS)
restored_ftot = restore_ftot(splits,vecs)
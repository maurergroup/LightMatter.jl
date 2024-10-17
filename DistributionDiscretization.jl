using Random,Distributions,JLD2,DelimitedFiles,OrdinaryDiffEq,DataInterpolations,Integrals,Roots
include("SimulationBuilder.jl")

function setup()
    las=define_laser_system(:Gaussian,fwhm=150,fluence=124.8,photon_en=1.55)
    sim = define_simulation_settings(nlelecphon=true,nlelecheat=true,noneqelec=true
    ,elecphonint=true,elecelecint=true,electemp=true,phonontemp=true)
    mp = define_material_parameters(las,extcof=14.9,gamma=6.117e-7,debye=343,noatoms=85,plasma=13.4,thermalcond=0.0025,
    elecperatom=1,eleceffmass=1.01,dos="DOS/Cu_DOS.dat",secmomspecfun=29e-6,elecphon=6.24e-7,ballistic=0.0,cph=0.015,τf=24.9)
    cons=Constants(8.617e-5,0.6582)
    dim = define_sim_dimensions(Dimension=1,Lengths=400,spacing=2)#define_sim_dimensions(Dimension=0)#
    return sim,mp,las,dim,cons
end

function generate_DOS(File::String,n)
    TotalDOS::Matrix{Float64}=readdlm(File)
    return DataInterpolations.LinearInterpolation(TotalDOS[:,2].*n,TotalDOS[:,1],extrapolate=true)
end

function discretization(x::Vector{Float64}, n::Int, w::Vector{Float64})
    @assert all(0 .<= x .<= 1) "All elements of x must be between 0 and 1"
    @assert length(x) == length(w) "x and w must have the same length"
    
    m = length(x)
    binary_vectors = zeros(Int, m, n) 

    for i in 1:m
        for j in 1:n
            binary_vectors[i, j] = rand(Bernoulli(x[i]))
        end
    end

    original_weighted_sum = sum(x .* w)
    tolerance = original_weighted_sum*0.01

    for j in 1:n
        current_weighted_sum = sum(binary_vectors[:, j] .* w)
        diff = original_weighted_sum - current_weighted_sum

        while abs(diff) > tolerance
            if diff > 0
                for i in 1:m
                    if abs(diff) <= tolerance
                        break
                    elseif binary_vectors[i, j] == 0 && w[i] > 0
                        binary_vectors[i, j] = 1
                        diff -= w[i]
                    end
                end
            elseif diff < 0
                for i in m:-1:1 
                    if abs(diff) <= tolerance
                        break
                    elseif binary_vectors[i, j] == 1 && w[i] > 0
                        binary_vectors[i, j] = 0
                        diff += w[i]
                    end
                end
            end
        end
    end
    return binary_vectors
end

function restore_ftot(splits,n)
    restore=zeros(length(splits[:,1]))
    for i in eachindex(splits[:,1])
        restore[i]=sum(splits[i,:])
    end
    return restore./n
end

@load "TestDistributions.jld2" sol

sim,mp,las,dim,cons = setup()

egrid = collect(range(-3*1.55,3*1.55,step=0.02))
DOS_spl = generate_DOS("DOS/Cu_DOS.dat",59)
DOS = DOS_spl(egrid)
fneq = sol[4:end,451]
Tel=sol[1,451]
n=sol[2,451]
μ=find_chemicalpotential(n,Tel,DOS_spl,cons.kB,mp.FE,mp.n0)
feq = FermiDirac(Tel,μ,cons.kB,egrid)
ftot = feq.+fneq

vecs=1000
println("Started")
splits = split_into_binary_vectors(ftot,vecs,DOS)
restored_ftot = restore_ftot(splits,vecs)
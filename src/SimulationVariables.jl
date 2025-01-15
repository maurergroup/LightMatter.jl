"""
    Determines the Fermi energy of the system. This is defined as the distance between the bottom and
    top of the valence band. This is used within Fermi Liquid Theory Relaxation Time where it is scaled
    based on the Fermi Energy^2. In all other places, the Fermi Energy is defined as 0.0.
"""
function get_FermiEnergy(File::String)
    TotalDOS::Matrix{Float64}=readdlm(File)
    Nonzero = findfirst(!=(0.0),TotalDOS[:,2])
    return abs(TotalDOS[Nonzero,1])
end
"""
    Converts a file location for the DOS into an interpolation object. It assumes that the DOS file
    is in units of states/atom and therefore scales the number of states by the number of atoms/nm(n).
"""
function generate_DOS(File::String,n)
    TotalDOS::Matrix{Float64}=readdlm(File)
    return get_interpolate(TotalDOS[:,1],TotalDOS[:,2].*n)
end

function spatial_DOS(folder::String,geometry::String,bulk::String,n::Real,dim::Dimension,tolerance)
    bulkDOS = readdlm(bulk)
    bulkDOSspl = Interpolations.interpolate(bulkDOS[:,1],bulkDOS[:,2]*n,SteffenMonotonicInterpolation())
    bulkDOSspl=extrapolate(bulkDOSspl,Flat())
    files,heights = get_files_heights_forDOS(folder,geometry,tolerance)
    DOS_1 = readdlm(folder*files[1],skipstart=4)
    egrid = DOS_1[:,1]
    zDOS = build_zDOSArray(egrid,folder,files,heights,bulkDOSspl,n)
    Temp=zeros(dim.length,length(egrid))
    for z in eachindex(dim.grid)
        for E in eachindex(egrid)
            Temp[z,E] = zDOS[E](dim.grid[z])
        end
    end
    #DOSScale(Temp,bulkDOSspl(egrid),egrid)
    zgridDOS=Vector{Any}(undef,dim.length)
    for i in eachindex(zgridDOS) 
        zgridDOS[i]=get_interpolate(egrid, Temp[i,:])
    end
    return zgridDOS
end

function DOSScale(Temp,bulk,Energies)
    fd = FermiDirac(0.0,0.0,8.617e-5,Energies)
    for i in eachindex(Temp[:,1])
        f(u) = extended_Bode(Energies,fd.*(u*Temp[i,:] .- bulk))
        x0 = 1
        prob = ZeroProblem(f,x0)
        rescale = solve(prob,Order16();atol=1e-10,rtol=1e-10)
        Temp[i,:] = Temp[i,:]*rescale
    end
    return Temp
end

function build_zDOSArray(egrid,folder,files,heights,bulkDOS,n)
    zDOS=Matrix{Float64}(undef,length(heights),length(egrid))
    for i in eachindex(files)
        TotalDOS=readdlm(folder*files[i],skipstart=4)
        zDOS[i,:]=TotalDOS[:,2]*n
    end
    zDOS=vcat(zDOS,transpose(bulkDOS(egrid)))
    zDOS=vcat(zDOS,transpose(bulkDOS(egrid)))
    heights=vcat(heights,findmax(heights)[1]+0.1)
    heights=vcat(heights,findmax(heights)[1]+0.1)
    zDOSspl=Vector{Any}(undef,length(egrid))
    for x in eachindex(zDOSspl)
        zDOSspl[x]=get_interpolate(heights,zDOS[:,x])
    end
    return zDOSspl
end

function get_files_heights_forDOS(folder,geometry,tolerance)
    files_from_folder = readdir(folder)
    dos_files = filter(f -> endswith(f, ".dat"), files_from_folder)
    split = splitext.(dos_files)
    file_names,extensions = [getindex.(split, i) for i in eachindex(first(split))]
    atoms = get_slabgeometry(geometry)
    layers = get_atomiclayers(atoms,tolerance)
    idxs = lpad.(convert.(Int,layers[:,1]),4,"0")
    files = []
    heights = []
    for i in eachindex(file_names)
        for j in eachindex(idxs)
            if endswith(file_names[i],idxs[j])
                push!(files,file_names[i]*extensions[j])
                push!(heights,layers[j,4])
            end
        end
    end
    heights = (heights .- heights[1])./10
    return files, heights
end

function get_atomiclayers(atoms,tolerance)
    unique_layers=[]
    for i in eachindex(atoms[:,1])
        to_push = true
        for j in eachindex(unique_layers)
            if abs(atoms[i,end]-unique_layers[j][end]) <= tolerance
                to_push = false
                break
            end
        end
        if to_push
            push!(unique_layers,atoms[i,:])
        end
    end
    return stack(unique_layers,dims=1)
end

function get_slabgeometry(file_path)
    atom_data = []
    top_constraint = [Inf,0]
    i=1
    geom = readdlm(file_path)
    for l in eachindex(geom[:,1])
        if geom[l,1] == "atom"
            if l != size(geom,1)
                if geom[l+1,1] != "constrain_relaxation"
                    push!(atom_data, [i, geom[l,2], geom[l,3], geom[l,4]])
                    i+=1
                end
            else
                push!(atom_data, [i, geom[l,2], geom[l,3], geom[l,4]])
                i+=1
            end
        end

    end
    return stack(atom_data,dims=1)
end
"""
    Generates an interpolation object with extrapolation from any two vectors of reals. The
    extrapolation returns the value of the boundaries. This should be suitable for DOS that are
    constant at the calculated boundaries and electronic distributions whose energy range is wide
    enough to capture all thermal and non-thermal behaviour.
"""
get_interpolate(xvals,yvals) = DataInterpolations.LinearInterpolation(yvals,xvals,extrapolate=true)
"""
    Sets up and solves the non-linear problem of determing the chemical potential at the current 
    electronic temperature.
"""
function find_chemicalpotential(no_part::Float64,Tel::Float64,DOS::spl,kB::Float64,egrid)::Float64
    f(u) = no_part - get_thermalparticles(u,Tel,DOS,kB,egrid)
    return solve(ZeroProblem(f,0.0),Order1();atol=1e-3,rtol=1e-3)
end

function get_thermalparticles(μ::Float64,Tel::Float64,DOS::spl,kB::Float64,egrid)::Float64
    return extended_Bode(DOS(egrid).*FermiDirac(Tel,μ,kB,egrid),egrid)
end
"""
    Determines the number of particles in any system using an interpolation of the system and
    the DOS of the system.
"""
function get_noparticles(Dis::Vector{<:Real},DOS::spl,egrid)
    return extended_Bode(Dis.*DOS(egrid),egrid)
end

function p_T(μ::Float64,Tel::Float64,DOS::spl,kB::Float64,egrid)
    return extended_Bode(dFDdT(kB,Tel,μ,egrid).*DOS(egrid),egrid)
end

function p_μ(μ::Float64,Tel::Float64,DOS::spl,kB::Float64,egrid)
    return extended_Bode(dFDdμ(kB,Tel,μ,egrid).*DOS(egrid),egrid)
end
"""
    Determines the internal energy of any system using an interpolation of that system and the
    DOS of the system.
"""
function get_internalenergy(Dis::Vector{<:Real},DOS::spl,egrid)
    return extended_Bode(Dis.*DOS(egrid).*egrid,egrid)
end

function c_T(μ::Float64,Tel::Float64,DOS::spl,kB::Float64,egrid::Vector{<:Real})
    return extended_Bode(dFDdT(kB,Tel,μ,egrid).*DOS(egrid).*egrid,egrid)
end

function c_μ(μ::Float64,Tel::Float64,DOS::spl,kB::Float64,egrid)
    return extended_Bode(dFDdμ(kB,Tel,μ,egrid).*DOS(egrid).*egrid,egrid)
end

function extended_Bode(y::Vector{<:Real}, x::Vector{<:Real})
    n = length(x)
    h = x[2]-x[1]  # The spacing between points
    integral = 0.0

    # Determine how many intervals to use for each rule
    integral_limit = n - 1

    # Apply Composite 3/8 Rule for the remaining intervals
    integral += (2h / 45) * (
        7 * y[1] + 7 * y[integral_limit+1] +
        32 * sum(y[2:2:integral_limit]) +
        12 * sum(y[3:4:integral_limit]) +
        14 * sum(y[4:3:integral_limit])
    )
    return integral
end
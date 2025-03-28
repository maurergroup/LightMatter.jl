"""
    Determines the Fermi energy of the system. This is defined as the distance between the bottom and
    top of the valence band. This is used within Fermi Liquid Theory Relaxation Time where it is scaled
    based on the Fermi Energy^2. In all other places, the Fermi Energy is defined as 0.0.
"""
function get_FermiEnergy(File)
    TotalDOS::Matrix{Float64}=readdlm(File)
    Nonzero = findfirst(!=(0.0),TotalDOS[:,2])
    return abs(TotalDOS[Nonzero,1])
end
"""
    Converts a file location for the DOS into an interpolation object. It assumes that the DOS file
    is in units of states/atom and therefore scales the number of states by the number of atoms/nm(n).
"""
function generate_DOS(File::String,unit_scalar)
    TotalDOS::Matrix{Float64}=readdlm(File)
    return get_interpolant(TotalDOS[:,1],TotalDOS[:,2]*unit_scalar)
end
"""
    Gets the unit cell volume as the bulk DOS in FHI-AIMS is in the units of states/eV/V_UC
"""
function get_unitcellvolume(geometry_file::String)
    geometry = readdlm(geometry_file)
    vectors = geometry[geometry[:,1] .== "lattice_vector",:] #Assumes FHI-aims geometry file
    a = vectors[1,2:4]
    b = vectors[2,2:4]
    c = vectors[3,2:4]
    return abs(dot(a,cross(b,c)))/1000 # converts Å^3 to nm^3
end
"""
    Overarching scheme for reading in a folder of atom projected DOS and returning a vector which matches that of
    Dimension.grid where each point in the vector is the DOS for that layer. This means taking the DOS from the folder
    interpolating in the Z-direction and extrapolating to bulk where necessary. Everything is currently 1D
"""
function spatial_DOS(folder::String,geometry::String,bulk::String,Vbulk,dim::Dimension,tolerance)
    bulkDOS = readdlm(bulk) #reads in the bulk DOS
    bulkDOSspl = get_interpolant(bulkDOS[:,1],bulkDOS[:,2]./Vbulk) #creates a spline for the bulk DOS
    files,heights = get_files_heights_forDOS(folder,geometry,tolerance) #get a vector of file names and their respective heights
    DOS_1 = readdlm(folder*files[1]) #Reads in a trial DOS 
    egrid = DOS_1[:,1]#Pulls the energy grid from the trial DOS as all folder DOS should be solved on same x-axis
    zDOS = build_zDOSArray(egrid,folder,files,heights)#Builds a vector in energy of splines of the DOS in the z direction
    Temp=zeros(dim.length,length(egrid)) #Temporary file to be filled with values from the interpolation vector above
    for z in eachindex(dim.grid)
        for E in eachindex(egrid)
            Temp[z,E] = zDOS[E](dim.grid[z]) #Calculates from the splines the values of each z and E point
        end
    end
    DOSScale(Temp,bulkDOSspl(egrid),egrid) #Scales all DOS' to the bulk dos to ensure particle conservation
    zgridDOS=Vector{Any}(undef,dim.length)
    for i in eachindex(zgridDOS) 
        zgridDOS[i]=get_interpolant(egrid, Temp[i,:]) # Builds the final array in the z-direction of splines of the DOS in energy
    end
    return zgridDOS
end
"""
    Ensures each DOS in the z-direction has the same number of particles as the bulk at 0K
"""
function DOSScale(Temp,bulk,Energies)
    fd = FermiDirac(0.0,0.0,Energies)
    for i in eachindex(Temp[:,1])
        f(u) = extended_Bode(u*fd.*Temp[i,:],Energies) - extended_Bode(fd.*bulk,Energies)
        x0 = 1
        prob = ZeroProblem(f,x0)
        rescale = solve(prob,Order16();atol=1e-12,rtol=1e-12)
        Temp[i,:] = Temp[i,:]*rescale
    end
    return Temp
end
"""
    Builds an array in the energy axis which contains a spline of the values in the z-axis
"""
function build_zDOSArray(egrid,folder,files,heights)
    zDOS=Matrix{Float64}(undef,length(heights),length(egrid))
    for i in eachindex(files)
        TotalDOS=readdlm(folder*files[i])
        zDOS[i,:]=TotalDOS[:,2]
    end
    zDOSspl=Vector{spl}(undef,length(egrid))
    for x in eachindex(zDOSspl)
        zDOSspl[x]=get_interpolant(heights,zDOS[:,x])
    end
    return zDOSspl
end
"""
    Returns a vector of the file names of each atom and their respective heights - assumes it's in Å 
"""
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
    heights = (heights .- heights[1])./10 #Å to nm
    return files, heights
end
"""
    Reduces the number of atoms to a single column where the tolerance decides the minimum distance
    for atoms to be deemed on seperate layers(defaults to 0.1 Å)
"""
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
"""
    Extracts the atoms and their positions from a geometry.in file
"""
function get_slabgeometry(file_path)
    atom_data = []
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
    get_interpolant(xvals::Vector{<:Real},yvals::Vector{<:Real})
    Generates an interpolation object with extrapolation from any two vectors of reals. The
    extrapolation returns the value of the boundaries. This should be suitable for DOS that are
    constant at the calculated boundaries and electronic distributions whose energy range is wide
    enough to capture all thermal and non-thermal behaviour. It is mainly used in the construction of DOS'
    but can be used elsewhere in the code when an interpolation is required.
"""
get_interpolant(xvals,yvals) = DataInterpolations.LinearInterpolation(yvals,xvals,extrapolation = ExtrapolationType.Constant)

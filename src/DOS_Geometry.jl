"""
    get_FermiEnergy(File::String)
    
    Calculates the Fermi energy defined as the difference between 0.0 and the
    bottom of the valence band. Assumes the DOS provided has the Fermi energy
    set to 0.0.

    # Arguments
    - 'File': Path to a file containing the density of states data.

    # Returns
    - The Fermi energy calculated at the bottom of the valence band.
"""
function get_FermiEnergy(File::String)
    TotalDOS = readdlm(File, comments=true)
    Nonzero = findfirst(!=(0.0), TotalDOS[:,2])
    return abs(TotalDOS[Nonzero,1])
end
"""
    generate_DOS(File::String, unit_scalar::Real)
    
    Generates a spline of a DOS from a file. Assumes the structure of the DOS is column 1 is Energy in eV
    and column 2 is States in eV⁻¹V⁻¹ (volume of unit cell) 

    # Arguments
    - 'File': Path to the total DOS file.
    - 'unit_scalar': Scalar to convert the units (1/V in nm⁻³).

    # Returns
    - An interpolation object representing the DOS.
"""
function generate_DOS(File::String, unit_scalar::Real)
    TotalDOS = readdlm(File,comments=true)
    return get_interpolant(TotalDOS[:,1], TotalDOS[:,2] * unit_scalar)
end
"""
    get_unitcellvolume(geometry_file::String)
    
    Calculates the volume of the unit cell for DOS unit conversion

    # Arguments
    - 'geometry_file': Path to the geometry.in file

    # Returns
    - Volume of the unit cell in nm⁻³
"""
function get_unitcellvolume(geometry_file::String)
    geometry = readdlm(geometry_file, comments=true)
    #atoms = count(x->x=="atom", geometry[:,1])
    vectors = geometry[geometry[:,1] .== "lattice_vector",:] #Assumes FHI-aims geometry file
    a = vectors[1,2:4]
    b = vectors[2,2:4]
    c = vectors[3,2:4]
    return (abs(dot(a,cross(b,c)))/1000) # converts Å^3 to nm^3
end
"""
    spatial_DOS(folder::String,geometry::String,bulk::String,Vbulk::Real,dim::Dimension,tolerance::Real)
    
    Creates a spline of a DOS at each z-grid point in the simulation. Reads a folder of atom projected DOS's and the 
    respective geomwtry.in file to determine the height of each DOS and interpolates between them to create the final
    DOS'

    # Arguments
    - 'folder': The folder where the atom-projected DOS' are present
    - 'geometry': The geometry.in file used in the creation of the atom-projected DOS'
    - 'bulk': The bulk DOS file location
    - 'Vbulk': The volume of the bulk unit cell
    - 'dim': The Dimension struct which holds the z-grid
    - 'tolerance': The minimum height 2 atoms need to be apart to be considered seperate layers (default = 0.1 nm)

    # Returns
    - A vector of splines corresponding to the DOS at each z-height
"""
function spatial_DOS(folder::String, geometry::String, bulk::String, Vbulk::Real, dim::Dimension, tolerance::Real)
    bulkDOS = readdlm(bulk, comments=true) #reads in the bulk DOS
    bulkDOSspl = get_interpolant(bulkDOS[:,1], bulkDOS[:,2] ./ Vbulk) #creates a spline for the bulk DOS
    files,heights = get_files_heights_forDOS(folder,geometry,tolerance) #get a vector of file names and their respective heights
    DOS_1 = readdlm(folder*files[1], comments=true) #Reads in a trial DOS 
    egrid = DOS_1[:,1] #Pulls the energy grid from the trial DOS as all folder DOS should be solved on same energy-axis
    zDOS = build_zDOSArray(egrid,folder, files, heights) #Builds a vector in energy of splines of the DOS in the z direction
    Temp = zeros(dim.length, length(egrid)) #Temporary file to be filled with values from the interpolation vector above
    for z in eachindex(dim.grid)
        for E in eachindex(egrid)
            Temp[z,E] = zDOS[E](dim.grid[z]) #Calculates from the splines the values of each z and E point
        end
    end
    DOSScale!(Temp,bulkDOSspl(egrid), egrid) #Scales all DOS' to the bulk dos to ensure particle conservation
    zgridDOS = Vector{spl}(undef, dim.length)
    for i in eachindex(zgridDOS) 
        zgridDOS[i] = get_interpolant(egrid, Temp[i,:]) # Builds the final array in the z-direction of splines of the DOS in energy
    end
    return zgridDOS
end
"""
    get_files_heights_forDOS(folder::String,geometry::String,tolerance::Real)
    
    Extracts atom from geometry, removes all bar one from each layer defined by tolerance, then connects the atom
    to the corresponding file in the folder of DOS'. Readjusts the heights to set the top layer to 0.0 and the
    lower layers to increase from there. 

    # Arguments
    - 'folder': The folder where the atom-projected DOS' are present
    - 'geometry': The geometry.in file used in the creation of the atom-projected DOS'
    - 'tolerance': The minimum height 2 atoms need to be apart to be considered seperate layers (default = 0.1Å)

    # Returns
    - A vector of DOS files and their respective heights
"""
function get_files_heights_forDOS(folder::String, geometry::String, tolerance::Real)
    files_from_folder = readdir(folder) # Reads all file names in folder
    dos_files = filter(f -> endswith(f, ".dat"), files_from_folder) #Filters out those that end .dat
    split = splitext.(dos_files) # Splits into matrix of file name ; extension
    file_names,extensions = [getindex.(split, i) for i in eachindex(first(split))] #Reformats split1
    atoms = get_slabgeometry(geometry) #Gets matrix of all atomic information, number and coordinate
    layers = get_atomiclayers(atoms, tolerance) #Removes all atoms other than 1 from each layer 
    idxs = lpad.(convert.(Int,layers[:,1]),4,"0") #Pads values to same length
    files = Vector{String}(undef, length(file_names))
    heights = zeros(length(file_names))
    for i in eachindex(file_names)
        for j in eachindex(idxs)
            if endswith(file_names[i],idxs[j])
                files[i] = file_names[i] * extensions[j]
                heights[i] = layers[j,4]
            end
        end
    end
    heights = (heights .- heights[1])./10 #Å to nm and sets the surface to 0.0
    return files, heights
end
"""
    get_slabgeometry(file_path::String)
    
    Extracts the atoms and their coordinates from a FHI-aims slab geometry.in file. It ignores any 
    atoms that have their relaxation constrained. 

    # Arguments
    - 'file_path': The geometry.in file used in the creation of the atom-projected DOS'

    # Returns
    - A matrix of the atom number and their coordinates
"""
function get_slabgeometry(file_path::String)
    atom_data = []
    i = 1
    geom = readdlm(file_path, comments=true)
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
    return stack(atom_data, dims=1)
end
"""
    get_atomiclayers(atoms::Matrix{<:Real},tolerance::Real)
    
    Seperats the atoms into their layers and selects a single atom from each layer. To remove degeneracy for 
    larger supercell structures.

    # Arguments
    - 'atoms': Matrix of the atom number and it's repseictve coordinates
    - 'tolerance': The minimum height 2 atoms need to be apart to be considered seperate layers (default = 0.1Å)

    # Returns
    - A trimmed matrix of atoms now containing one atom per layer
"""
function get_atomiclayers(atoms::Matrix{<:Real}, tolerance::Real)
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
    return stack(unique_layers, dims=1)
end
"""
    build_zDOSArray(egrid::Vector{<:Real},folder::String,files::Vector{String},heights::Vector{<:Real})
    
    Builds a matrix of the DOS as a function of height and energy for the individual layers. 

    # Arguments
    - 'egrid': Energy grid the DOS is calculated on
    - 'folder': The folder where the atom-projected DOS' are present
    - 'files': Vector of file names 
    - 'heights': Vector of each file names height

    # Returns
    - A matrix of states as a function of height and energy
"""
function build_zDOSArray(egrid::Vector{<:Real}, folder::String, files::Vector{String}, heights::Vector{<:Real})
    zDOS = Matrix{Float64}(undef, length(heights), length(egrid))
    for i in eachindex(files)
        TotalDOS = readdlm(folder*files[i], comments=true)
        zDOS[i,:] = TotalDOS[:,2]
    end
    zDOSspl = Vector{spl}(undef, length(egrid))
    for x in eachindex(zDOSspl)
        zDOSspl[x] = get_interpolant(heights, zDOS[:,x])
    end
    return zDOSspl
end
"""
    DOSScale!(Temp::Matrix{<:Real},bulk::Vector{<:Real},Energies::Vector{<:Real})
    
    Ensures that all DOS are scaled to the same number of particles as the bulk

    # Arguments
    - 'Temp': Matrix of number of states in height X Energy
    - 'bulk': The bulk DOS on the same energy grid as Temp
    - 'Energies': The energy grid used for Temp and bulk

    # Returns
    - Temp with the rescaled DOS'
"""
function DOSScale!(Temp::Matrix{<:Real}, bulk::Vector{<:Real}, Energies::Vector{<:Real})
    fd = FermiDirac(0.0,0.0, Energies)
    for i in eachindex(Temp[:,1])
        f(u) = extended_Bode(u*fd.*Temp[i,:], Energies) - extended_Bode(fd.*bulk, Energies)
        x0 = 1
        prob = ZeroProblem(f,x0)
        rescale = solve(prob, Order16(); atol=1e-12, rtol=1e-12)
        Temp[i,:] = Temp[i,:] * rescale
    end
    return Temp
end
"""
    get_interpolant(xvals::Vector{<:Real},yvals::Vector{<:Real})
    
    Generates a linear spline of any two vectors with a constant extrapolation applied to the boundaries.

    # Arguments
    - 'xvals': x-axis of the desired spline
    - 'yvals': y-axis of the desired spline

    # Returns
    - Spline of yvals vs xvals
"""
@inline get_interpolant(xvals::Vector{<:Real}, yvals::Vector{<:Real}) = DataInterpolations.LinearInterpolation(yvals, xvals, extrapolation = ExtrapolationType.Constant)
"""
    DOS_initialization(bulk_DOS::Union{String,Vector{String}}, bulk_geometry::String, DOS_folder::String, slab_geometry::String,
                       atomic_layer_tolerance::Real, dimension::Dimension, zDOS::Bool, DOS::Union{Nothing, spl})
    
    Determines the desired DOS configuration and assembles it accordingly

    # Arguments
    - 'bulk_DOS': The bulk DOS file location
    - 'bulk_geometry': The bulk DOS geometry location
    - 'DOS_folder': The folder where the atom-projected DOS' are present
    - 'slab_geometry': The geometry.in file used in the creation of the atom-projected DOS'
    - 'atomic_layer_tolerance': The minimum height 2 atoms need to be apart to be considered seperate layers (default = 0.1 nm)
    - 'dim': The Dimension struct which holds the z-grid
    - 'zDOS': Bool for enabling a spatially dependent DOS
    ' 'DOS': Allows the user to use their own splined DOS, for other regions of the code it must be the same type as normal DOS
    
    # Returns
    - A spline or vector of splines for the desired DOS structure
"""
function DOS_initialization(bulk_DOS::Union{String,Vector{String}}, bulk_geometry::String, DOS_folder::Union{Nothing,String}, slab_geometry::Union{Nothing,String},
                            atomic_layer_tolerance::Real, dimension::Dimension, zDOS::Bool, DOS::Union{Nothing, spl})
    if DOS !== nothing
        return DOS
    else
        if bulk_DOS isa String
            Vbulk = get_unitcellvolume(bulk_geometry)
            if zDOS == true 
                DOS = spatial_DOS(DOS_folder, slab_geometry, bulk_DOS, Vbulk, dimension, atomic_layer_tolerance)
            else
                DOS = generate_DOS(bulk_DOS, 1/Vbulk) 
            end
        else
            Vbulk = zeros(length(bulk_DOS))
            DOS = Vector{spl}(undef, length(bulk_DOS))
            for i in eachindex(bulk_DOS)
                Vbulk[i] = get_unitcellvolume(bulk_geometry[i])
                if zDOS == true 
                    DOS[i] = spatial_DOS(DOS_folder[i], slab_geometry[i], bulk_DOS[i], Vbulk[i],dimension, atomic_layer_tolerance)
                else
                    DOS[i] = generate_DOS(bulk_DOS[i], 1/Vbulk[i])
                end
            end
        end
        return DOS
    end
end
"""
    DOS_initialization(bulk_DOS::Union{String,Vector{String}}, bulk_geometry::String, DOS_folder::String, slab_geometry::String,
                       atomic_layer_tolerance::Float64, dimension::Dimension, zDOS::Bool, DOS::Union{Nothing, spl})
    
    Determines the desired DOS configuration and assembles it accordingly. Returns missing if no DOS files are provided.
    Spatial DOS with a μ offset for antenna-reactor complexes is currently not supported.

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
    - A spline or vector of splines for the desired DOS structure, or missing if no files are provided
"""
function DOS_initialization(bulk_DOS::Union{String,Vector{String},Nothing}, bulk_geometry::Union{String,Vector{String},Nothing},
                            DOS_folder::Union{Nothing,String}, slab_geometry::Union{Nothing,String},
                            atomic_layer_tolerance::Float64, dimension::Dimension, zDOS::Bool, DOS::Union{Nothing, spl})
    if DOS !== nothing
        return DOS
    else
        if bulk_DOS isa Nothing
            return missing
        elseif bulk_DOS isa String
            if bulk_geometry isa Nothing
                return missing
            end
            Vbulk = get_unitcellvolume(bulk_geometry)
            if zDOS == true 
                DOS = spatial_DOS(DOS_folder, slab_geometry, bulk_DOS, Vbulk, dimension, atomic_layer_tolerance)
            else
                DOS = generate_DOS(bulk_DOS, 1/Vbulk)
            end
        else
            Vbulk = zeros(length(bulk_DOS))
            DOS = Vector{spl}(undef, length(bulk_DOS))
            offset = calculate_μoffset(DOS, Vbulk)
            for i in eachindex(bulk_DOS)
                if bulk_geometry isa Nothing
                    return missing
                end
                Vbulk[i] = get_unitcellvolume(bulk_geometry[i])
                if zDOS == true 
                    DOS[i] = spatial_DOS(DOS_folder[i], slab_geometry[i], bulk_DOS[i], Vbulk[i], dimension, atomic_layer_tolerance)
                else
                    DOS[i] = generate_DOS(bulk_DOS[i], 1/Vbulk[i]; offset=offset[i])
                end
            end
        end
        return DOS
    end
end
"""
    get_FermiEnergy(File::String)
    
    Extracts the Fermi energy from a DOS file, finding the first non-zero DOS value

    # Arguments
    - 'File': Path to the total DOS file in eV format

    # Returns
    - The Fermi energy in eV
"""
function get_FermiEnergy(File::String)
    TotalDOS = readdlm(File, comments=true)
    Nonzero = findfirst(!=(0.0), TotalDOS[:,2])
    return abs(TotalDOS[Nonzero,1])
end
"""
    generate_DOS(File::String, unit_scalar::Float64)
    
    Generates a spline of a DOS from a file. Assumes the structure of the DOS is column 1 is Energy in eV
    and column 2 is States in eV⁻¹V⁻¹ (volume of unit cell) 

    # Arguments
    - 'File': Path to the total DOS file.
    - 'units_scalar': Scalar to convert the units (1/V in nm⁻³).

    # Returns
    - An interpolation object representing the DOS.
"""
function generate_DOS(File::String, units_scalar::Float64;offset=0.0)
    TotalDOS = readdlm(File,comments=true)
    return get_interpolant(TotalDOS[:,1].+offset, TotalDOS[:,2] * units_scalar)
end
"""
    build_DOS(dos_file::String, geometry_file::String)

    A convenient constructor for building a DOS from the DOS and geometry file.

    # Arguments
    - 'dos_file': Path to the total DOS file.
    - 'geometry_file': Path to the geometry.in file

    # Returns
    - An interpolation object representing the DOS.
"""
function build_DOS(dos_file::String, geometry_file::String)
    v = get_unitcellvolume(geometry_file::String)
    return generate_DOS(dos_file, 1/v)
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
    vectors = geometry[geometry[:,1] .== "lattice_vector",:] #Assumes FHI-aims geometry file
    a = vectors[1,2:4]
    b = vectors[2,2:4]
    c = vectors[3,2:4]
    return (abs(dot(a,cross(b,c)))/1000) # converts Å^3 to nm^3
end

function get_atomicdensity(geometry_file::String)
    geometry = readdlm(geometry_file, comments=true)
    atoms = count(x->contains(x, "atom"), geometry[:,1])
    vectors = geometry[geometry[:,1] .== "lattice_vector",:] #Assumes FHI-aims geometry file
    a = vectors[1,2:4]
    b = vectors[2,2:4]
    c = vectors[3,2:4]
    return atoms / (abs(dot(a,cross(b,c)))/1000) # converts Å^3 to nm^3
end
"""
    spatial_DOS(folder::String,geometry::String,bulk::String,Vbulk::Float64,dim::Dimension,tolerance::Float64)
    
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
function spatial_DOS(folder::String, geometry::String, bulk::String, Vbulk::Float64, dim::Dimension, tolerance::Float64)
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
    get_files_heights_forDOS(folder::String,geometry::String,tolerance::Float64)
    
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
function get_files_heights_forDOS(folder::String, geometry::String, tolerance::Float64)
    files_from_folder = readdir(folder) # Reads all file names in folder
    dos_files = filter(f -> endswith(f, ".dat"), files_from_folder) #Filters out those that end .dat
    split = splitext.(dos_files) # Splits into matrix of file name ; extension
    file_names,extensions = [getindex.(split, i) for i in eachindex(first(split))] #Reformats split1
    atoms = get_slabgeometry(geometry) #Gets matrix of all atomic information, number and coordinate
    layers = get_atomiclayers(atoms, tolerance) #Removes all atoms other than 1 from each layer 

    files = Vector{String}(undef, size(layers, 1))
    heights = zeros(size(layers, 1))
    for i in eachindex(layers[:,1])
        if endswith(file_names[i],"$i")
            files[i] = file_names[i] * extensions[i]
            heights[i] = layers[i,4]
        end
    end
    heights = (heights .- heights[1])./10 #Å to nm and sets the surface to 0.0
    perm = sortperm(heights)
    return files[perm], heights[perm]
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
                end
            else
                push!(atom_data, [i, geom[l,2], geom[l,3], geom[l,4]])
            end
        end

    end
    stk_data = stack(atom_data, dims=1)
    stk_data[:,4] .-= stk_data[1,4]
    return stk_data
end
"""
    get_atomiclayers(atoms::Matrix{Float64},tolerance::Float64)
    
    Seperats the atoms into their layers and selects a single atom from each layer. To remove degeneracy for 
    larger supercell structures.

    # Arguments
    - 'atoms': Matrix of the atom number and it's repseictve coordinates
    - 'tolerance': The minimum height 2 atoms need to be apart to be considered seperate layers (default = 0.1Å)

    # Returns
    - A trimmed matrix of atoms now containing one atom per layer
"""
function get_atomiclayers(atoms::Matrix{Float64}, tolerance::Float64)
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
    build_zDOSArray(egrid::Vector{Float64},folder::String,files::Vector{String},heights::Vector{Float64})
    
    Builds a matrix of the DOS as a function of height and energy for the individual layers. 

    # Arguments
    - 'egrid': Energy grid the DOS is calculated on
    - 'folder': The folder where the atom-projected DOS' are present
    - 'files': Vector of file names 
    - 'heights': Vector of each file names height

    # Returns
    - A matrix of states as a function of height and energy
"""
function build_zDOSArray(egrid::Vector{Float64}, folder::String, files::Vector{String}, heights::Vector{Float64})
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
    DOSScale!(Temp::Matrix{Float64},bulk::Vector{Float64},Energies::Vector{Float64})
    
    Ensures that all DOS are scaled to the same number of particles as the bulk

    # Arguments
    - 'Temp': Matrix of number of states in height X Energy
    - 'bulk': The bulk DOS on the same energy grid as Temp
    - 'Energies': The energy grid used for Temp and bulk

    # Returns
    - Temp with the rescaled DOS'
"""
function DOSScale!(Temp::Matrix{Float64}, bulk::Vector{Float64}, Energies::Vector{Float64})
    fd = FermiDirac(0.0,0.0, Energies)
    for i in eachindex(Temp[:,1])
        f(u,p) = integration_algorithm(u*fd.*Temp[i,:], Energies) - integration_algorithm(fd.*bulk, Energies)
        x0 = 1
        prob = NonlinearProblem(f,x0)
        rescale = solve(prob, SimpleNewtonRaphson(); atol=1e-12, rtol=1e-12)
        Temp[i,:] = Temp[i,:] * rescale.u
    end
    return Temp
end
"""
    get_interpolant(xvals::Vector{Float64},yvals::Vector{Float64})
    
    Generates a linear spline of any two vectors with a constant extrapolation applied to the boundaries.

    # Arguments
    - 'xvals': x-axis of the desired spline
    - 'yvals': y-axis of the desired spline

    # Returns
    - Spline of yvals vs xvals
"""


@inline get_interpolant(xvals, yvals) = Interpolations.linear_interpolation(xvals, yvals, extrapolation_bc = Flat())#Interpolations.interpolate((xvals,), yvals, Gridded(Linear()))
"""
    build_group_velocity(v_g::Union{Vector{Float64},Nothing}, FE::Union{Float64,Vector{Float64}}, Conductivity::Bool, conductive_velocity::Symbol, structure::Structure)
    
    Creates a vector or array of vectors (spatial DOS) for the group veolcity for ballistic electron transport. Users can also provide a constant value in the form
    of v_g, they must also set conductive_veolcity to constant.
    Returns missing if DOS is missing or if conductivity is disabled.
    Currently Implemented:
    - :fermigas : Assumes a free electron gas solution therefore is an analytical form of the group velocity
    - :effectiveoneband : Uses the effective one band model to convert a DOS into a group velocity, for more details see Mueller & Rethfeld, Phys. Rev. B 87, 035139.
    - :constant : Uses the v_g argument to set a constant group velocity for all energy ranges

    # Arguments
    - 'v_g': A constant group velocity value if :constant is requested
    - 'FE': The Fermi energy, calculated from get_FermiEnergy
    - 'Conductivity': Sets whether ballistic transport should be enabled
    - 'conductive_velocity': The form the user wants the group velocity to take
    - 'structure': Contains all structural information including DOS and number of elemental systems

    # Returns
    - The group velocity vector or array of vectors as requested by the user for ballistic electron transport, or missing if not applicable
"""
function build_group_velocity(v_g::Union{Vector{Float64},Missing,Float64}, FE::Union{Float64,Vector{Float64},Missing}, Conductivity::Bool, conductive_velocity::Symbol, structure::Structure{1})
    Conductivity || return missing

    if ismissing(v_g)
        return build_group_velocity_missing(structure, FE, conductive_velocity)
    elseif conductive_velocity == :constant
        return convert_units(u"nm/fs", v_g)
    end

    return missing
end

function build_group_velocity(v_g::Union{Vector{Float64},Missing,Float64}, FE::Union{Float64,Vector{Float64},Missing}, Conductivity::Bool, conductive_velocity::Symbol, structure::Structure{N}) where {N}
    Conductivity || return missing

    if ismissing(v_g)
        return build_group_velocity_missing(structure, FE, conductive_velocity)
    elseif conductive_velocity == :constant
        return build_group_velocity_constant(v_g, structure)
    end

    return missing
end

function build_group_velocity(v_g::Union{Vector{Float64},Missing,Float64}, FE::AbstractVector, Conductivity::Bool, conductive_velocity::Symbol, structure::Structure{N}) where {N}
    Conductivity || return missing

    if ismissing(v_g) || ismissing(structure.DOS)
        return missing
    elseif conductive_velocity == :constant
        return build_group_velocity_constant(v_g, structure)
    elseif conductive_velocity == :fermigas
        return build_group_velocity_fermigas(structure, FE)
    elseif conductive_velocity == :effectiveoneband
        return build_group_velocity_effective(structure, FE)
    end

    return missing
end

function build_group_velocity_missing(structure::Structure{1}, FE, conductive_velocity::Symbol)
    if ismissing(structure.DOS) || ismissing(FE)
        return missing
    elseif conductive_velocity == :fermigas
        return get_fermigas_velocity(Ref(structure.egrid), FE)
    elseif conductive_velocity == :effectiveoneband
        if structure.Spatial_DOS == true
            vg = zeros(structure.dimension.length, length(structure.egrid))
            for i in 1:structure.dimension.length
                vg[i, :] = effective_one_band_velocity(structure.bandstructure[i].k_to_E, structure.DOS[i], structure.egrid, FE)
            end
            return vg
        else
            return effective_one_band_velocity(structure.bandstructure.k_to_E, structure.DOS, structure.egrid, FE)
        end
    end

    return missing
end

function build_group_velocity_missing(structure::Structure{N}, FE::AbstractVector, conductive_velocity::Symbol) where {N}
    if ismissing(structure.DOS)
        return missing
    elseif conductive_velocity == :fermigas
        return build_group_velocity_fermigas(structure, FE)
    elseif conductive_velocity == :effectiveoneband
        return build_group_velocity_effective(structure, FE)
    end

    return missing
end

function build_group_velocity_constant(v_g, structure::Structure{1})
    return convert_units(u"nm/fs", v_g)
end

function build_group_velocity_constant(v_g, structure::Structure{N}) where {N}
    vg = Vector{Vector{Float64}}(undef, N)
    grids = split_grid(structure.dimension.grid, structure.dimension.InterfaceHeight)
    for i in eachindex(grids)
        vg[i] = fill(convert_units(u"nm/fs", v_g), length(grids[i]))
    end
    return vg
end

function build_group_velocity_fermigas(structure::Structure{N}, FE) where {N}
    vg = Vector{Vector{Float64}}(undef, N)
    grids = split_grid(structure.dimension.grid, structure.dimension.InterfaceHeight)
    for i in eachindex(grids)
        vg[i] = fill(get_fermigas_velocity(Ref(structure.egrid), FE[i]), length(grids[i]))
    end
    return vg
end

function build_group_velocity_effective(structure::Structure{N}, FE) where {N}
    if structure.Spatial_DOS == true
        vg = zeros(structure.dimension.length, length(structure.egrid))
        for i in 1:structure.dimension.length
            vg[i, :] = effective_one_band_velocity(structure.bandstructure[i].k_to_E, structure.DOS[i], structure.egrid, FE[i])
        end
        return vg
    else
        vg = zeros(structure.dimension.length, length(structure.egrid))
        Threads.@threads for i in eachindex(vg[:, 1])
            j = mat_picker(structure.dimension.grid[i], structure.dimension.InterfaceHeight)
            vg[i, :] .= effective_one_band_velocity(structure.bandstructure[j].k_to_E, structure.DOS[j], structure.egrid, FE[i])
        end
        return vg
    end
end
"""
    get_fermigas_velocity(egrid::Vector{Float64}, EF::Float64)
    
    The analytical free electron gas group velocity, requested by conductive_velocity = :fermigas

    # Arguments
    - 'egrid': Energy grid all distributions are solved on
    - 'FE': The Fermi energy, calculated from get_FermiEnergy

    # Returns
    - The free electron gas group velocity
"""
function get_fermigas_velocity(egrid::Vector{Float64}, EF::Float64)
    return sqrt.(2 * (egrid.+EF) ./ Constants.me)
end
"""
    get_fermigas_dos(egrid, FE)
    
    Function for calculating a free electron gas.

    # Arguments
    - 'egrid': Energy grid all distributions are solved on
    - 'FE': The Fermi energy, calculated from get_FermiEnergy

    # Returns
    - The free electron gas denisty-of-states
"""
function get_fermigas_dos(egrid::Vector{Float64}, EF::Float64)
    comp1 = 1/(3*π^2)
    comp2 = (2*Constants.me/Constants.ħ^2)^(3/2)
    comp3 = (egrid.+FE).^0.5
    return comp1 .* comp2 .* comp3
end
"""
    effective_one_band_velocity(DOS::spl, egrid::Vector{Float64}, FE::Float64)
    
    Calculates the group velocity from the effective one band model.
    For more details see Mueller & Rethfeld, Phys. Rev. B 87, 035139.

    # Arguments
    - 'DOS': The density-of-states of the system
    - 'egrid': Energy grid all distributions are solved on
    - 'FE': The Fermi energy, calculated from get_FermiEnergy

    # Returns
    - The effective one band model group velocity as a vector or vector of vectors depending on the structure
    of the system
"""
function effective_one_band_velocity(Bandstructure::Spline1D, DOS::spl, egrid::Vector{Float64}, FE::Float64)
    k_E = effective_onebandmodel(DOS,egrid,FE)
    v_g = similar(k_E)
    if eltype(v_g) <: AbstractVector
        for j in eachindex(v_g)
            dE_dk = Dierckx.derivative.(Ref(Bandstructure),k_E)
            kE_spl = LightMatter.get_interpolant(egrid,k_E)
            dEdk_spl = LightMatter.get_interpolant(k_E,dE_dk)
            v_g[j] = dEdk_spl(kE_spl(egrid))./ Constants.ħ
        end
    else
        dE_dk = Dierckx.derivative.(Ref(Bandstructure),k_E)
        kE_spl = LightMatter.get_interpolant(egrid,k_E)
        dEdk_spl = LightMatter.get_interpolant(k_E,dE_dk)
        v_g = dEdk_spl(kE_spl(egrid))./ Constants.ħ
    end
    return v_g
end
"""
    effective_onebandmodel(DOS::spl, egrid::Vector{Float64}, FE::Float64)
    
    Calculates the dispersion relation within the effective one band model.
    For more details see Mueller & Rethfeld, Phys. Rev. B 87, 035139.

    # Arguments
    - 'DOS': The density-of-states of the system
    - 'egrid': Energy grid all distributions are solved on
    - 'FE': The Fermi energy, calculated from get_FermiEnergy

    # Returns
    - The effective one band model dispersion relation
"""
function effective_onebandmodel(DOS, egrid::Vector{Float64}, FE::Float64)
    k_E = zeros(length(egrid))
    
    factor = 6π
    int(u,p) = DOS(u)

    for (i,E) in enumerate(egrid)
        prob=IntegralProblem(int, -FE, E)
        sol = solve(prob, HCubatureJL(initdiv=100), abstol=1e-8, reltol=1e-8)
        k_E[i] = cbrt(factor*sol.u)
    end
    return k_E
end

"""
    bandstructure_initialization(bandstructure, DOS, egrid, FE)
    
    Initializes the band+
     structure interpolation objects for the effective one band model

    # Arguments
    - 'bandstructure': Type of band structure model (:effectiveoneband or other)
    - 'DOS': Density of states (can be spl or Vector of spl)
    - 'egrid': Energy grid
    - 'FE': Fermi energy (can be Float64 or Vector for multi-element systems)

    # Returns
    - Spline objects for band structure interpolation, or missing if not applicable
"""
function bandstructure_initialization(Bandstructure, DOS, egrid, FE)
    if Bandstructure == nothing
        if ismissing(DOS)
            return missing
        end
        if !(typeof(DOS) <: spl)
            E_k = Vector{bandstructure}(undef, length(DOS))
            for i in eachindex(DOS)
                if length(DOS) == length(FE)
                    fe = FE[i]
                else
                    fe = FE
                end
                temp_k = effective_onebandmodel(DOS[i], egrid, fe)
                E_k[i] = bandstructure(Dierckx.Spline1D(temp_k, egrid, k=3, bc="nearest"),
                                      Dierckx.Spline1D(egrid, temp_k, k=3, bc="nearest"))
            end
            return E_k
        else
            temp_k = effective_onebandmodel(DOS, egrid, FE)
            return bandstructure(Dierckx.Spline1D(temp_k, egrid, k=3, bc="nearest"),
                                 Dierckx.Spline1D(egrid, temp_k, k=3, bc="nearest"))
        end
    else
        return Bandstructure
    end
end

"""
    dos_velocity(DOS::spl, egrid::Vector{Float64}, FE::Float64)
    
    Calculates the group velocity directly from the density-of-states without assuming 
    a parabolic bandstructure. Uses the effective one-band model inversion:
    
    v_g(E) = (ℏ/3) * (3π * ∫ g(E') dE')^(2/3) / g(E)
    
    This avoids reconstructing a parabolic E(k) and gives more realistic velocities
    for materials with non-parabolic bands (e.g., gold).
    For more details see Mueller & Rethfeld, Phys. Rev. B 87, 035139.

    # Arguments
    - 'DOS': The density-of-states spline
    - 'egrid': Energy grid all distributions are solved on
    - 'FE': The Fermi energy

    # Returns
    - Group velocity v_g(E) as a vector in units of nm/fs (if ℏ is eV⋅fs)
"""
function dos_velocity(DOS::spl, egrid::Vector{Float64}, FE::Float64)
    # Calculate the integrated DOS at each energy
    factor = 6π
    int(u, p) = DOS(u)
    
    integrated_dos = zeros(length(egrid))
    for (i, E) in enumerate(egrid)
        prob = IntegralProblem(int, -FE, E)
        sol = solve(prob, HCubatureJL(initdiv=100), abstol=1e-8, reltol=1e-8)
        integrated_dos[i] = sol.u
    end
    
    # Calculate v_g(E) = (ℏ/3) * (3π * N(E))^(2/3) / g(E)
    dos_values = DOS.(egrid)
    v_g = similar(egrid)
    
    for i in eachindex(egrid)
        if dos_values[i] > 1e-10  # Avoid division by zero
            v_g[i] = (Constants.ħ / 3) * (factor * integrated_dos[i])^(2/3) / dos_values[i]
        else
            v_g[i] = 0.0
        end
    end
    
    return v_g
end

function calculate_μoffset(DOS::Vector{String}, Vbulk)
    grid = collect(range(-20.0,20.0, length=1001))
    tmp_DOS = Vector{Spl}(undef, length(DOS))
    N0s = zeros(length(DOS))
    Dmax = zeros(length(DOS))
    offset = zeros(length(DOS))
    for i in eachindex(DOS)
        ref_line = findall(x->x[13] == "mu", eachrow(readdlm(DOS[i])))
        offset[i] = reference[ref_line[1],15]
        tmp_DOS[i] = generate_DOS(DOS[i], 1/Vbulk[i], offset=offset[i])
        Dmax[i] = integration_algorithm(tmp_DOS[i], grid)
        N0s[i] = get_thermalparticles(0.0, 1e-16, tmp_DOS[i], grid)
    end
    Ntotal = sum(N0s)

    function root_function(μ)
        X = sum(get_thermalparticles(μ, 1e-16, tmp_DOS[i], grid) for i in eachindex(tmp_DOS))
        return X - Ntotal
    end

    Emin = grid[1]
    Emax = grid[end]

    margin = 1e-6

    lo = Emin - margin
    hi = Emax + margin

    flo = root_function(lo)
    fhi = root_function(hi)

    # Expand bracket if needed.
    niter = 0
    while flo > 0 && niter < 100
        margin *= 2
        lo -= margin
        flo = root_function(lo)
        niter += 1
    end

    niter = 0
    while fhi < 0 && niter < 100
        margin *= 2
        hi += margin
        fhi = root_function(hi)
        niter += 1
    end

    if sign(flo) == sign(fhi)
        error("""
        Could not bracket common chemical potential.

        root_function(lo) = $flo at lo = $lo
        root_function(hi) = $fhi at hi = $hi

        Check DOS energy ranges and normalization.
        """)
    end

    mu_common = bisect_root(root_function, lo, hi;)

    return offset .- mu_common   
end

function bisect_root(f, a, b; tol = 1e-10, maxiter = 200)
    fa = f(a)
    fb = f(b)

    if abs(fa) < tol
        return a
    elseif abs(fb) < tol
        return b
    end

    if sign(fa) == sign(fb)
        error("Root is not bracketed: f(a) = $fa, f(b) = $fb")
    end

    for _ in 1:maxiter
        c = 0.5 * (a + b)
        fc = f(c)

        if abs(fc) < tol || abs(b - a) < tol
            return c
        end

        if sign(fc) == sign(fa)
            a = c
            fa = fc
        else
            b = c
            fb = fc
        end
    end

    return 0.5 * (a + b)
end
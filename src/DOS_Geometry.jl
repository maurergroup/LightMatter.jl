"""
    DOS_initialization(bulk_DOS::Union{String,Vector{String}}, bulk_geometry::String, DOS_folder::String, slab_geometry::String,
                       atomic_layer_tolerance::Float64, dimension::Dimension, zDOS::Bool, DOS::Union{Nothing, spl})
    
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
function DOS_initialization(bulk_DOS::Union{String,Vector{String},Nothing}, bulk_geometry::Union{String,Vector{String},Nothing},
                            DOS_folder::Union{Nothing,String}, slab_geometry::Union{Nothing,String},
                            atomic_layer_tolerance::Float64, dimension::Dimension, zDOS::Bool, DOS::Union{Nothing, spl})
    if DOS !== nothing
        return DOS
    else
        if bulk_DOS isa Nothing
            return get_interpolant([1,2,3], [4,5,6])
        elseif bulk_DOS isa String
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
    generate_DOS(File::String, unit_scalar::Float64)
    
    Generates a spline of a DOS from a file. Assumes the structure of the DOS is column 1 is Energy in eV
    and column 2 is States in eV⁻¹V⁻¹ (volume of unit cell) 

    # Arguments
    - 'File': Path to the total DOS file.
    - 'unit_scalar': Scalar to convert the units (1/V in nm⁻³).

    # Returns
    - An interpolation object representing the DOS.
"""
function generate_DOS(File::String, unit_scalar::Float64)
    TotalDOS = readdlm(File,comments=true)
    return get_interpolant(TotalDOS[:,1], TotalDOS[:,2] * unit_scalar)
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
    atoms = count(x->x=="atom", geometry[:,1])
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
        f(u,p) = Bode_rule(u*fd.*Temp[i,:], Energies) - Bode_rule(fd.*bulk, Energies)
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
@inline get_interpolant(xvals, yvals) = Interpolations.linear_interpolation(xvals, yvals, extrapolation_bc = Flat())
"""
    build_group_velocity(v_g::Union{Vector{Float64},Nothing}, FE::Union{Float64,Vector{Float64}}, Conductivity::Bool, conductive_velocity::Symbol, structure::Structure)
    
    Creates a vector or array of vectors (spatial DOS) for the group veolcity for ballistic electron transport. Users can also provide a constant value in the form
    of v_g, they must also set conductive_veolcity to constant.
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
    - The group velocity vector or array of vectors as requested by the user for ballistic electron transport
"""
function build_group_velocity(v_g::Union{Vector{Float64},Nothing}, FE::Union{Float64,Vector{Float64}}, Conductivity::Bool, conductive_velocity::Symbol, structure::Structure)
    if isnothing(v_g)
        if Conductivity == true
            if conductive_velocity == :fermigas
                if structure.Elemental_System > 1
                    elements = structure.Elemental_System
                    v_g = Vector{Vector{Float64}}(undef,elements)
                    grids = split_grid(structure.dimension.grid,structure.dimension.InterfaceHeight)
                    for i in 1:Elemental_System
                        length = length(grids[i])
                        v_g[i] = fill(get_fermigas_velocity(Ref(structure.egrid),FE),length)
                    end
                else
                    return get_fermigas_velocity(Ref(structure.egrid),FE)
                end
            elseif conductive_velocity == :effectiveoneband
                if structure.Spatial_DOS == true
                    v_g = zeros(structure.dimension.length,length(structure.egrid))
                    for i in 1:structure.dimension.length
                        v_g[i,:] = effective_one_band_velocity(structure.bandstructure[i][1], structure.DOS[i],structure.egrid,FE[i])
                    end
                elseif structure.Elemental_System != 1
                    v_g = zeros(structure.dimension.length,length(structure.egrid))
                    Threads.@threads for i in eachindex(v_g[:,1])
                        j = mat_picker(structure.dimension.grid[i], structure.dimension.InterfaceHeight)
                        v_g[i,:] .= effective_one_band_velocity(structure.bandstructure[j][1], structure.DOS[j],structure.egrid,FE[j])
                    end
                    return v_g
                else
                    v_g = effective_one_band_velocity(structure.bandstructure[1], structure.DOS,structure.egrid,FE)
                end
            end
        else 
            return [NaN]
        end
    else
        if Conductivity == true
            if conductive_velocity == :constant
                if structure.Elemental_System > 1
                    elements = structure.Elemental_System
                    V_g = Vector{Vector{Float64}}(undef,elements)
                    grids = split_grid(structure.dimension.grid,structure.dimension.InterfaceHeight)
                    for i in 1:Elemental_System
                        length = length(grids[i])
                        V_g[i] = fill(convert_units(u"nm/fs", v_g),length)
                    end
                else
                    return convert_units(u"nm/fs", v_g)
                end
            end
        else 
            return [NaN]
        end
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
function effective_one_band_velocity(bandstructure::AkimaInterpolation, DOS::spl, egrid::Vector{Float64}, FE::Float64)
    k_E = effective_onebandmodel(DOS,egrid,FE)
    v_g = similar(k_E)
    if eltype(v_g) <: AbstractVector
        for j in eachindex(v_g)
            dE_dk = DataInterpolations.derivative.(Ref(bandstructure),k_E)
            kE_spl = LightMatter.get_interpolant(egrid,k_E)
            dEdk_spl = LightMatter.get_interpolant(k_E,dE_dk)
            v_g[j] = dEdk_spl(kE_spl(egrid))./ Constants.ħ
        end
    else
        dE_dk = DataInterpolations.derivative.(Ref(bandstructure),k_E)
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
    
    factor = 3π^2
    int(u,p) = DOS(u)

    for (i,E) in enumerate(egrid)
        prob=IntegralProblem(int, -FE, E)
        sol = solve(prob, HCubatureJL(initdiv=100), abstol=1e-8, reltol=1e-8)
        k_E[i] = cbrt(factor*sol.u)
    end
    return k_E
end

function bandstructure_initialization(bandstructure, DOS, egrid, FE)
    if bandstructure == :effectiveoneband
        if !(typeof(DOS) <: spl)
            E_k = Vector{Vector{AkimaInterpolation}}(undef, length(DOS))
            for i in eachindex(DOS)
                if length(DOS) == length(FE)
                    fe = FE[i]
                else
                    fe = FE
                end
                temp_k = effective_onebandmodel(DOS[i], egrid, fe)

                E_k[i] = [DataInterpolations.AkimaInterpolation(egrid,temp_k,extrapolation = ExtrapolationType.Constant),
                          DataInterpolations.AkimaInterpolation(temp_k, egrid,extrapolation = ExtrapolationType.Constant)]
            end
            return E_k
        else
            temp_k = effective_onebandmodel(DOS, egrid, FE)
            return [DataInterpolations.AkimaInterpolation(egrid,temp_k,extrapolation = ExtrapolationType.Constant),
                    DataInterpolations.AkimaInterpolation(temp_k, egrid,extrapolation = ExtrapolationType.Constant)]
        end
    else
        return [DataInterpolations.AkimaInterpolation([1,2,3],[4,5,6],extrapolation = ExtrapolationType.Constant),
                DataInterpolations.AkimaInterpolation([1,2,3],[4,5,6],extrapolation = ExtrapolationType.Constant)]
    end
end
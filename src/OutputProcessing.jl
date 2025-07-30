###
# Add keywords to save additional information - change the output::Symbol system
# Open and close the file regularly so things don't break so easily
# This generally needs a good tidying
###
"""
    post_production(sol,file_name::String,initial_temps::Dict{String,Float64},output::Symbol,sim::Simulation)
    
    Handles processing, saving the simulation after it has completed. Uses HDF5 file format
    Currently the only output setting supported is :minimum which only saves the parameters, chemical potential and
    progated systems

    # Arguments
    - 'sol': The solution of the simulation
    - 'file_name': The name the user wishes the file to be saved as 
    - 'initial_temps': Initial temperature the baths have been set to
    - 'output': The output method required - may change in future
    - 'sim': Simulation settings and parameters

    # Returns
    - Nothing is returned but a file is created
"""
function post_production(sol, file_name::String, initial_temps::Dict{String,Float64}, output::Vector{Symbol}, sim::Simulation)
    temp_name = "temp_"*file_name[1:end-5]*".jld2" 
    @save temp_name sol
    fid = create_datafile_and_structure(file_name)
    write_simulation(fid,sim::Simulation)

    results = seperate_results(sol,initial_temps,sim)
    fid["Time Span"] = results["times"]
    write_dynamicalvariables(fid, results)

    selected_output_functions(fid, results, sim, output)

    create_group(fid, "System Equations")
    output_functions(fid["System Equations"], sim)
    rm(temp_name)
    close(fid)
end
"""
    create_datafile_and_structure(file_name::String)
    
    Creates a file with the standard data structure:
    - Athermal Electrons
    - Density Matrix
    - Electronic Temperature
    - Phononic Temperature
    - Electronic Distribution
    - Phononic Distribution
    - Laser
    - Structure

    # Arguments
    - 'file_name': The name the user wishes the file to be saved as 

    # Returns
    - A created HDF5 file
"""
function create_datafile_and_structure(file_name::String)
    fid = h5open(file_name, "w")
    create_group(fid, "Athermal Electrons")
    create_group(fid, "Density Matrix")
    create_group(fid, "Electronic Temperature")
    create_group(fid, "Phononic Temperature")
    create_group(fid, "Electronic Distribution")
    create_group(fid, "Phononic Distribution")
    create_group(fid, "Laser")
    create_group(fid, "Structure")
    return fid
end
"""
    dict_to_hdf5(f,d::Dict{String,Any})
    
    Unpacks a dictionary and saves it to the designated file location (f)

    # Arguments
    - 'f': Location in the HDF5 file the dictionary information will be saved to
    - 'd': The dictionary to be unpacked and saved

    # Returns
    - Nothing returned but dictionary saved to location
"""
function dict_to_hdf5(f, d)
    for (key, value) in d
        if typeof(value) <:Vector{<:Vector}
            value = stack(value, dims=1)
        end
        f[key] = value
    end
end
"""
    convert_symbols_to_strings(dict::Dict{Any,Any})
    
    Converts any Symbol key values to their respecitve string for saving

    # Arguments
    - 'dict': The dictionary which needs keys converting

    # Returns
    - Dictionary with all Symbol keys changed to strings
"""
function convert_symbols_to_strings(dict)
    return Dict(
        k => (v isa Symbol ? String(v) : v) 
        for (k, v) in dict
    )
end
"""
    write_simulation(f,sim::Simulation)
    
    Writes all settings and parameters to their respective location in the file (f)

    # Arguments
    - 'f': File to be written to
    - 'sim': Simulation settings and parameters

    # Returns
    - Nothing returned but all parameters and settings wrote to file
"""
function write_simulation(f,sim::Simulation)
    athermal = Dict{String,Any}(String(key)=>getfield(sim.athermalelectrons, key) for key ∈ fieldnames(AthermalElectrons))
    dict_to_hdf5(f["Athermal Electrons"], convert_symbols_to_strings(athermal))

    electronic_t = Dict{String,Any}(String(key)=>getfield(sim.electronictemperature, key) for key ∈ fieldnames(ElectronicTemperature))
    dict_to_hdf5(f["Electronic Temperature"], convert_symbols_to_strings(electronic_t))

    phononic_t = Dict{String,Any}(String(key)=>getfield(sim.phononictemperature, key) for key ∈ fieldnames(PhononicTemperature))
    dict_to_hdf5(f["Phononic Temperature"], convert_symbols_to_strings(phononic_t))

    electronic_d = Dict{String,Any}(String(key)=>getfield(sim.electronicdistribution, key) for key ∈ fieldnames(ElectronicDistribution))
    dict_to_hdf5(f["Electronic Distribution"], convert_symbols_to_strings(electronic_d))

    phononic_d = Dict{String,Any}(String(key)=>getfield(sim.phononicdistribution, key) for key ∈ fieldnames(PhononicDistribution))
    dict_to_hdf5(f["Phononic Distribution"], convert_symbols_to_strings(phononic_d))

    density_m = Dict{String,Any}(String(key)=>getfield(sim.densitymatrix, key) for key ∈ fieldnames(DensityMatrix))
    dict_to_hdf5(f["Density Matrix"], convert_symbols_to_strings(density_m))

    laser = Dict{String,Any}(String(key)=>getfield(sim.laser, key) for key ∈ fieldnames(Laser))
    dict_to_hdf5(f["Laser"], convert_symbols_to_strings(laser))
    
    extract_structure(f, sim.structure)
end
"""
    extract_structure(f, structure::Structure)
    
    Specific way of writing structure settings to the file as both DOS and Dimension need special handling

    # Arguments
    - 'f': File to be written to
    - 'structure': Structure settings and parameters

    # Returns
    - Nothing returned but structure parameters and settings wrote to file
"""
function extract_structure(f, structure::Structure)
    create_group(f["Structure"], "Dimension")
    struc = Dict("Elemental_System" => structure.Elemental_System, "Spatial_DOS" => structure.Spatial_DOS, "egrid" => structure.egrid)
    merge!(struc, write_DOS(structure))

    dimension = Dict{String,Any}(String(key) => getfield(structure.dimension, key) for key ∈ fieldnames(Dimension))
    dict_to_hdf5(f["Structure"]["Dimension"], convert_symbols_to_strings(dimension))
    dict_to_hdf5(f["Structure"], convert_symbols_to_strings(struc))
end
"""
    write_DOS(structure::Structure)
    
    Handles converting the DOS splines into a writable format

    # Arguments
    - 'structure': Structure settings and parameters

    # Returns
    - Dictionary of the DOS
"""
function write_DOS(structure::Structure)
    
    egrid = collect(range(-20,20, step=0.01))
    if typeof(structure.DOS) == Vector{spl} && structure.Elemental_System == 1
        DOS = zeros(length(structure.DOS), length(egrid))
        for i in eachindex(structure.DOS)
            DOS[i,:] = structure.DOS[i](egrid)
        end
    elseif structure.Elemental_System > 1 && structure.DOS isa Vector
        DOS = zeros(length(structure.DOS), length(egrid))
        for i in eachindex(structure.DOS)
            DOS[i,:] = structure.DOS[i](egrid)
        end
    else
        DOS = zeros(length(egrid))
        DOS[:] = structure.DOS(egrid)
    end
    return Dict("DOS" => DOS)
end
"""
    seperate_results(sol, initial_temps::Dict{String,Float64}, sim::Simulation)
    
    Seperates the results held inside of the solution object
    Also fills all unpropagated subsystems with parameter/temp information where neccessary

    # Arguments
    - 'sol': The solution of the simulation
    - 'initial_temps': Initial temperature the baths have been set to
    - 'sim': Simulation settings and parameters

    # Returns
    - Dictionary of values of each of the seperated systems 
"""
function seperate_results(sol, initial_temps::Dict{String,Float64}, sim::Simulation)
    fields = propertynames(sol[1])
    vals = generate_valuedict(sol, sim, fields)
    populate_value_dict!(sol, fields, vals)
    populate_unpropagatedvalues!(initial_temps, fields, sim, vals)
    merge!(vals,Dict("times" => sol.t))
    return vals
end
"""
    generate_valuedict(sol, sim::Simulation, fields::Vector{Symbol})
    
    Generates a dictionary of all propagated systems with placeholder arrays

    # Arguments
    - 'sol': The solution of the simulation
    - 'sim': Simulation settings and parameters
    - 'fields': Symbolic name of the propgated subsystems

    # Returns
    - Dictionary of values of each of the propgated subsystems 
"""
function generate_valuedict(sol, sim::Simulation, fields)
    vals = Dict{String,Union{Float64,AbstractArray}}()
    for i in fields
        if i in [:Tel, :Tph, :noe]
            merge!(vals,Dict(String(i) => zeros(length(sol.t), sim.structure.dimension.length)))
        elseif i == :fneq
            merge!(vals,Dict(String(i) => zeros(length(sol.t), sim.structure.dimension.length, length(sim.structure.egrid))))
        end
    end
    return vals
end
"""
    populate_value_dict!(sol ,fields::Vector{Symbol}, vals::Dict{String,AbstractArray{Float64}})
    
    Populates the subsystem dictionary with the reuslting vlaues from the simulation

    # Arguments
    - 'sol': The solution of the simulation
    - 'fields': Symbolic name of the propgated subsystems
    - 'vals': Dictionary of the propgated subsytems with placeholder arrays

    # Returns
    - The vals dictionary with the actual results inside
"""
function populate_value_dict!(sol, fields, vals)
    for t in eachindex(sol.t)
        array = ArrayPartition(sol[t])
        for f in eachindex(fields)
            vals[String(fields[f])][t,:,:] .= array.x[f]
        end
    end
    return vals
end
"""
    populate_unpropagatedvalues!(sol, initial_temps::Dict{String,Float64}, fields::Vector{Symbol}, sim::Simulation, vals::Dict{String,AbstractArray{Float64}})
    
    Adds placeholder information to any unpropagated fields

    # Arguments
    - 'sol': The solution of the simulation
    - 'initial_temps': Initial temperature the baths have been set to
    - 'sim': Simulation settings and parameters
    - 'fields': Symbolic name of the propgated subsystems
    - 'vals': Dictionary of the propgated subsytems with placeholder arrays

    # Returns
    - vals dictionary with the added unpropagated subsystems
"""
function populate_unpropagatedvalues!(initial_temps::Dict{String,Float64}, fields, sim::Simulation, vals)
    if :Tel ∉ fields
        merge!(vals, Dict("Tel" => fill(initial_temps["Tel"], sim.structure.dimension.length)'))
    end
    if :Tph ∉ fields
        merge!(vals, Dict("Tph" => fill(initial_temps["Tph"], sim.structure.dimension.length)'))
    end
    if :noe ∉ fields
        merge!(vals, Dict("noe" => particlenumber_values(sim)))
    end
    return vals
end
"""
    particlenumber_values(sim::Simulation)
    
    Calculates the number of thermal electrons for a system when not propagated

    # Arguments
    - 'sim': Simulation settings and parameters

    # Returns
    - Array of the number of thermal electrons at each spatial point
"""
function particlenumber_values(sim::Simulation)
    L = sim.structure.dimension.length
    if typeof(sim.structure.DOS) == Vector{spl} && length(sim.structure.DOS) == L #DOS per z point
        no_part = zeros(L)
        for j in eachindex(sim.structure.dimension.grid)
            no_part[j] = get_thermalparticles(0.0, 1e-32, sim.structure.DOS[j], sim.structure.egrid)
        end
    elseif typeof(sim.structure.DOS) == Vector{spl} #DOS for each material rather than heights
        no_part = zeros(L)
        for j in eachindex(sim.structure.dimension.grid)
            mat = mat_picker(sim.structure.dimension.grid[j], sim.structure.dimension.InterfaceHeight)
            no_part[j] = get_thermalparticles(0.0,1e-32, sim.structure.DOS[mat], sim.structure.egrid)
        end
    else #One DOS
        no_part = fill(get_thermalparticles(0.0, 1e-32, sim.structure.DOS,sim.structure.egrid), L)
    end
    return no_part'
end

function write_dynamicalvariables(f, results)
    for (name, result) in results
        if name == "Tel"
            write_dataset(f["Electronic Temperature"],"Temperature",result)
        elseif  name == "Tph"
            write_dataset(f["Phononic Temperature"],"Temperature",result)
        elseif name == "fneq"
            write_dataset(f["Athermal Electrons"],"Non-Equilibrium Distribution",result)
        elseif name =="noe"
            write_dataset(f,"Particle Number",result)
        end
    end
end

function selected_output_functions(f, results, sim, output)
    for sym in output
        func_name = Symbol("output_", sym)
        if isdefined(LightMatter, func_name)
            func = getfield(LightMatter, func_name)
            func(f, results, sim,)
        else
            @warn "Function $func_name not defined"
        end
    end
end

function output_chemicalpotential(f, results, sim)
    if !haskey(results, "μ")
        Tel = results["Tel"]
        cp = zeros(size(Tel))
        dos = sim.structure.DOS
        n = results["noe"]

        Threads.@threads for i in eachindex([:,1])
            DOS = get_DOS(dos, sim, i)
            cp[i,:] .= find_chemicalpotential.(n[i,:],Tel[i,:],Ref(DOS),Ref(sim.structure.egrid))
        end

        write_dataset(f, "Chemical Potential", cp)
        merge!(results, Dict("μ" => cp))
    end
end

function output_ThermalFermiDistribution(f, results, sim)
    Tel = results["Tel"]
    FD=zeros(length(Tel[:,1]), sim.structure.dimension.length, length(sim.structure.egrid))
    output_chemicalpotential(f, results, sim)
    μ = results["μ"]
    Threads.@threads for i in eachindex(Tel[:,1])
        for j in eachindex(sim.structure.dimension.grid)
            FD[i,j,:] .= FermiDirac(Tel[i,j], μ[i,j], sim.structure.egrid)
        end
    end
    write_dataset(f["Electronic Temperature"], "Electronic Distribution", FD)
    merge!(results, Dict("TFD" => FD))
end 

function output_dTeldt(f, results, sim)
    Tel = results["Tel"]
    dTemp = similar(Tel)
    Threads.@threads for i in eachindex(dTemp[:,1])
        if i ==1
            dTemp[i,:] .= zeros(length(dTemp[i,:]))
        else
            dTemp[i,:] .= Tel[i,:] .- Tel[i-1,:]
        end
    end
    write_dataset(f["Electronic Temperature"], "Temperature Derivative", dTemp)
end

function output_TelConductivity(f, results, sim)
    Tel = results["Tel"]
    Tph = results["Tph"]
    spat = similar(Tel)
    Threads.@threads for i in eachindex(Tel[:,1])
        @views electrontemperature_conductivity!(Tel[i,:], sim.electronictemperature.κ, sim.structure.dimension.spacing, Tph[i,:], spat[i,:])
    end
    write_dataset(f["Electronic Temperature"], "Thermal Conductivity", spat)
end

function output_electronphononcoupling(f, results, sim)
    Tel = results["Tel"]
    Tph = results["Tph"]
    eps = similar(Tel)
    if sim.electronictemperature.ElectronPhononCouplingValue == :variable
        (; ω, λ) = sim.electronictemperature
        output_chemicalpotential(f, results, sim)
        μ = results["μ"]
        dos = sim.structure.DOS
        Threads.@threads for i in eachindex(Tel[1,:])
            DOS = get_DOS(dos, sim, i)
            (omega, lambda) = get_parameterhvalue((ω, λ), sim, i)
            @views eps[i,:] .= variable_electronphononcoupling.(omega, lambda, Ref(DOS), Tel[:,i], μ[:,i], Tph[:,i])
        end
    else
        Threads.@threads for i in eachindex(Tel[1,:])
            g_val = get_parameterhvalue(sim.electronictemperature.g, sim, i)
            @views eps[i,:] .= g_val * (Tel[:,i].-Tph[:,i])
        end
    end
    write_dataset(f["Electronic Temperature"], "Electron-Phonon Coupling", eps)
end     

function output_electronheatcapacity(f, results, sim)
    Tel = results["Tel"]
    hc = similar(Tel)
    if sim.electronictemperature.ElectronicHeatCapacity == nonlinear
        output_chemicalpotential(f, results, sim)
        μ = results["μ"]
        dos = sim.structure.DOS
        Threads.@threads for i in eachindex(Tel[1,:])
            DOS = get_DOS(dos, sim, i)
            @views hc[i,:] .= nonlinear_electronheatcapacity.(Tel[:,i], μ[:,i], Ref(DOS))
        end
    else
        Threads.@threads for i in eachindex(Tel[:,1])
            gamma =  get_parameterhvalue(sim.electronictemperature.γ, sim, i)
            @views hc[i,:] .= gamma*Tel[i,:]
        end
    end
    write_dataset(f["Electronic Temperature"], "Electronic Heat Capacity", hc)
end

function output_athermalelectronelectronscattering(f, results, sim)
    if !haskey(results, "e*erelax")
        fneq = results["fneq"]
        tel = results["Tel"]
        n = results["noe"]
        output_chemicalpotential(f, results, sim)
        μ = results["μ"]
        rel = similar(fneq)
        τ_expr = athem_relaxationtime(sim)
        lifetime = mk_function((sim,μ,Tel,),(),τ_expr)
        grid = sim.structure.dimension.grid
        sims = vec_simulation(sim)

        Threads.@threads for i in eachindex(fneq[:,1,1])
            for j in eachindex(tel[1,:])
                X = mat_picker(grid[j], sim.structure.dimension.InterfaceHeight)
                τee = lifetime(sims[X], μ[i,j], tel[i,j])
                @views rel[1,j,:] .= -athem_electronelectronscattering(tel[i,j], μ[i,j], sim, fneq[i,j,:], DOS, n[i,j], τee)
            end
        end
        write_dataset(f["Athermal Electrons"], "Athermal Electron-Electron Scattering", rel)
        merge!(results, Dict("e*erelax" => rel))
    end
end

function output_athermalelectronelectronenergychange(f, results, sim)
    output_athermalelectronelectronscattering(f, results, sim)
    relax = results["e*erelax"]
    uee = similar(results["Tel"])
    dos = sim.structure.DOS
    Threads.@threads for i in eachindex(uee[:,1])
        for j in eachindex(uee[i,:])
            DOS = get_DOS(dos, sim, j)
            uee[i,j] = get_internalenergy(-1*relax[:,j,i],Ref(DOS), Ref(sim.structure.egrid))
        end
    end
    write_dataset(f["Electronic Temperature"], "Athermal Electron-Electron Energy Flow", uee)
end

function output_phononicheatcapacity(f, results, sim)
    if sim.phononictemperature.PhononicHeatCapacity == :variable
        tph = results["Tph"]
        cph = similar(tph)
        (; n, θ) = sim.phononictemperature
        Threads.@threads for i in eachindex(Tph[:,1])
            for j in eachindex(Tph[1,:])
                (no_elec, θD) = get_parameterhvalue((n, θ), sim, j)
                cph[i,j] .= variable_phononheatcapacity(tph[i,j], no_elec, θD)
            end
        end
        write_dataset(f["Phononic Temperature"], "Phonon Heat Capacity", cph)
    else
        write_dataset(f["Phononic Temperature"], "Phonon Heat Capacity", sim.phononictemperature.Cph)
    end
end

function output_athemelectronphononscattering(f, results, sim)
    if !haskey(results, "e*phrelax")
        fneq = results["fneq"]
        Tel = results["Tel"]
        Tph = results["Tph"]
        tmp = similar(fneq)
        τ_expr = phonon_relaxationtime(sim)
        lifetime = mk_function((sim,Tel,Tph),(),τ_expr)
        grid = sim.structure.dimension.grid
        sims = vec_simulation(sim)

        Threads.@threads for i in eachindex(Tel[:,1])
            for j in eachindex(uep[1,:])
                X = mat_picker(grid[j], sim.structure.dimension.InterfaceHeight)
                @views @. tmp[i,j,:] = -fneq[i,j,:] / lifetime(sims[X], Tel[i,j], Tph[i,j])
            end
        end
        write_dataset(f["Athermal Electrons"], "Athermal Electron-Phonon Scattering", tmp)
        merge!(results, Dict("e*phrelax" => tmp))
    end
end

function output_athermalelectronphononenergychange(f, results, sim)
    output_athermalelectronphononscattering(f, results, sim)
    relax = results["e*phrelax"]
    uep = similar(results["Tel"])
    dos = sim.structure.DOS
    Threads.@threads for i in eachindex(uep[:,1])
        for j in eachindex(uep[i,:])
            DOS = get_DOS(dos, sim, j)
            uep[i,j] = get_internalenergy(-1*relax[:,j,i],Ref(DOS), Ref(sim.structure.egrid))
        end
    end
    write_dataset(f["Phononic Temperature"], "Athermal Electron-Phonon Energy Flow", uee)
end

function output_dndt(f, results, sim)
    n = results["n"]
    δn = similar(n)
    δn[1,:] = zeros(length(δn[1,:]))
    Threads.@threads for i in 2:length(n[:,1])
        @views @. δn[i,:] = n[i,:] .- n[i-1,:]
    end
    write_dataset(f["Electronic Temeprature"], "Change in Thermal Particle Number", δn)
end

function output_dfneqdt(f, results, sim)
    fneq = results["fneq"]
    dfneq = zeros(size(fneq))
    dfneq[1,:,:] .= zeros(size(fneq[1,:,:]))
    Threads.@threads for i in 2:length(fneq[:,1,1])
        @views @. dfneq = fneq[i,:,:] - fneq[i-1,:,:]
    end
    write_dataset(f["Athermal Electrons"], "Change in fneq", dfneq)
end

function output_athemexcitation(f, results, sim) #Currently working on this and below
    fneq = results["fneq"]
    Tel = results["Tel"]
    output_chemicalpotential(f, results, sim)
    μ = results["μ"]
    dos = sim.structure.DOS
    grid = sim.structure.dimension.grid
    sims = vec_simulation(sim)

    excite = similar(fneq)
    M = mk_function((), (), athem_excitation_matrixelements(sim))
    las = mk_function((t, sim, i),(), laser_factory(sim))
    Threads.@threads for i in eachindex(fneq[:,1,1])
        for j in eachindex(fneq[1,:,1])
            X = mat_picker(grid[j], sim.structure.dimension.InterfaceHeight)
            ftot = fneq[i,j,:] .+ FermiDirac(Tel[i,j], μ[i,j], sim.structure.egrid)
            DOS = get_DOS(dos, sim, j)
            @views excite[i,j,:] .= las(t, sims[X], j) * athemexcitation(ftot, sim.structure.egrid, DOS, sim.laser.hv, M())
        end
    end
    write_dataset(f["Athermal Electrons"], "Excitation", excite)
end 

function output_laserprofile(f, results, sim)
    time = results["times"]
    grid = sim.structure.dimension.grid
    temp_prof = zeros(length(time),length(grid))
    las = mk_function((t, sim, i),(), laser_factory(sim))

    sims = vec_simulation(sim)

    for i in eachindex(time)
        for j in eachindex(grid)
            X = mat_picker(grid[j], sim.structure.dimension.InterfaceHeight)
            temp_prof[i,j] = las(time[i], sims[X], j)
        end
    end
    write_dataset(f["Laser"], "Temporal Profile", temp_prof)
end

function write_dataset(file,dataset,data)
    dims = size(data)

    if ndims(data) == 1
        chunk_size = size(data)  # Safe chunking for 1D
    elseif ndims(data) == 2
        chunk_size = min.(dims, (1, dims[2]))  # Safe chunking for 2D
    elseif ndims(data) == 3
        chunk_size = min.(dims, (1, dims[2], dims[3]))  # Safe chunking for 3D
    else
        error("Unsupported data dimensionality: $(ndims(data))D")
    end

    HDF5.write_dataset(file, dataset, data, chunk=chunk_size, shuffle=true, deflate=3)
end

function get_DOS(DOS, sim, i)
    if DOS isa AbstractArray
        if length(DOS) == sim.structure.dimension.length
            return DOS[i]
        else 
            X = mat_picker(sim.structure.dimension.grid[i],sim.structure.dimension.InterfaceHeight)
            return DOS[X]
        end
    else
        return DOS
    end
end   

function get_parameterhvalue(val, sim, i)
    if sim.structure.Elemental_System > 1
        X = mat_picker(sim.structure.dimension.grid[i],sim.structure.dimension.InterfaceHeight)
        return getindex(val, X)
    else
        return val
    end
end

function get_parameterhvalue(tup::Tuple, sim, i)
    if sim.structure.Elemental_System > 1
        X = mat_picker(sim.structure.dimension.grid[i],sim.structure.dimension.InterfaceHeight)
        return getindex.(tup, X)
    else
        return tup
    end
end

function vec_simulation(sim::Simulation)
    if sim.structure.Elemental_System > 1
        sims = sim_seperation(sim)
    else 
        sims = sim
    end
    return sims
end

function output_functions(f,sim)
    sys = function_builder(sim)
    for i in keys(sys)
        f[i] = string(sys[i])
    end
    simulation_expr = simulation_construction(sys,sim)
    f["Total Equation Block"] = string(simulation_expr)
end
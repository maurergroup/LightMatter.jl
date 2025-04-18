###
# Add keywords to save additional information - change the output::Symbol system
# Open and close the file regularly so things don't break so easily
# This generally needs a good tidying
###
"""
    post_production(sol,file_name::String,initial_temps::Dict{String,<:Real},output::Symbol,sim::Simulation)
    
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
function post_production(sol, file_name::String, initial_temps::Dict{String,<:Real}, output::Symbol, sim::Simulation)
    fid = create_datafile_and_structure(file_name)
    write_simulation(fid,sim::Simulation)

    results = seperate_results(sol,initial_temps,sim)
    fid["Time Span"] = results["times"]
    μs = pp_chemicalpotential(results["Tel"], results["noe"], sim)
    fid["Chemical Potential"] = μs
    FD = pp_FermiDistribution(results["Tel"],sim,μs)


    create_group(fid, "System Equations")
    output_functions(fid["System Equations"], sim)

    if output==:Full
        println("Not Implemented running minimum")
        write_minimum(fid,results,FD,sim)
        close(fid)
    elseif output==:Minimum
        write_minimum(fid, results, FD, sim)
        close(fid)
    else 
        println("The only output options are full and minimum. Files cannot be easily overwritten so the file has been deleted")
        close(fid)
        rm(file_name)
    end
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
function dict_to_hdf5(f, d::Dict{String,Any})
    for (key, value) in d
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
function convert_symbols_to_strings(dict::Dict{Any,Any})
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
    elseif structure.Elemental_System > 1
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
    seperate_results(sol, initial_temps::Dict{String,<:Real}, sim::Simulation)
    
    Seperates the results held inside of the solution object
    Also fills all unpropagated subsystems with parameter/temp information where neccessary

    # Arguments
    - 'sol': The solution of the simulation
    - 'initial_temps': Initial temperature the baths have been set to
    - 'sim': Simulation settings and parameters

    # Returns
    - Dictionary of values of each of the seperated systems 
"""
function seperate_results(sol, initial_temps::Dict{String,<:Real}, sim::Simulation)
    fields = propertynames(sol[1])
    vals = generate_valuedict(sol, sim, fields)
    populate_value_dict!(sol, fields, vals)
    populate_unpropagatedvalues!(sol, initial_temps, fields, sim, vals)
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
function generate_valuedict(sol, sim::Simulation, fields::Vector{Symbol})
    vals = Dict{String,Union{Real,AbstractArray}}()
    for i in fields
        if i in [:Tel, :Tph, :noe]
            if i == :noe && sim.athermalelectrons.Enabled == false
                nothing
            else
                merge!(vals,Dict(String(i) => zeros(length(sol.t), sim.structure.dimension.length)))
            end
        elseif i == :fneq
            merge!(vals,Dict(String(i) => zeros(length(sol.t), sim.structure.dimension.length, length(sim.structure.egrid))))
        end
    end
    return vals
end
"""
    populate_value_dict!(sol ,fields::Vector{Symbol}, vals::Dict{String,AbstractArray{<:Real}})
    
    Populates the subsystem dictionary with the reuslting vlaues from the simulation

    # Arguments
    - 'sol': The solution of the simulation
    - 'fields': Symbolic name of the propgated subsystems
    - 'vals': Dictionary of the propgated subsytems with placeholder arrays

    # Returns
    - The vals dictionary with the actual results inside
"""
function populate_value_dict!(sol ,fields::Vector{Symbol}, vals::Dict{String,AbstractArray{<:Real}})
    for t in eachindex(sol.t)
        array = ArrayPartition(sol[t])
        for f in eachindex(fields)
            vals[String(fields[f])][t,:,:] .= array.x[f]
        end
    end
    return vals
end
"""
    remove_TTM_explicit_electrons!(vals::Dict{String,AbstractArray{<:Real}}, fields::Vector{Symbol}, sim::Simulation)
    
    WIP!
    Removes the explicit electron number from the value dictionary for the TTM region when the embedding method is used

    # Arguments
    - 'vals': Dictionary of the propgated subsytems with placeholder arrays
    - 'fields': Symbolic name of the propgated subsystems
    - 'sim': Simulation settings and parameters

    # Returns
    - vals dictionary with the thermal electron number array reduced in size
"""
function remove_TTM_explicit_electrons!(vals::Dict{String,AbstractArray{<:Real}}, fields::Vector{Symbol}, sim::Simulation)
    if sim.ParameterApprox.EmbeddingMethod == true
        for i in [:fneq,:noe]
            if i in fields
                vals[String(i)] = vals[String(i)][:,1,:]
            end
        end
    end
    return vals
end
"""
    populate_unpropagatedvalues!(sol, initial_temps::Dict{String,<:Real}, fields::Vector{Symbol}, sim::Simulation, vals::Dict{String,AbstractArray{<:Real}})
    
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
function populate_unpropagatedvalues!(sol, initial_temps::Dict{String,<:Real}, fields::Vector{Symbol}, sim::Simulation, vals::Dict{String,AbstractArray{<:Real}})
    if :Tel ∉ fields
        merge!(vals, Dict("Tel" => fill(initial_temps["Tel"], sim.structure.dimension.length)))
    end
    if :Tph ∉ fields
        merge!(vals, Dict("Tph" => fill(initial_temps["Tph"], sim.structure.dimension.length)))
    end
    if :noe ∉ fields
        merge!(vals, Dict("noe" => particlenumber_values(sim) ))
    end
    if :fneq ∉ fields
        merge!(vals, Dict("fneq" => zeros(length(sol.t), sim.structure.dimension.length, length(sim.structure.egrid))))
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
    if sim.athermalelectrons.Enabled == true && sim.athermalelectrons.AthermalElectron_ElectronCoupling == false
        if typeof(sim.structure.DOS) == Vector{spl} && sim.structure.Elemental_System == 1
            no_part = zeros(sim.structure.dimension.length)
            for j in eachindex(sim.structure.dimension.grid)
                no_part[j] = get_thermalparticles(0.0, 1e-32, sim.structure.DOS[j], sim.structure.egrid)
            end
        elseif typeof(sim.structure.DOS) == Vector{spl} && sim.structure.Elemental_System > 1
            no_part = zeros(sim.structure.dimension.length)
            for j in eachindex(sim.structure.dimension.grid)
                mat = mat_picker(sim.structure.dimension.grid[j], sim.structure.dimension.InterfaceHeight)
                no_part[j] = get_thermalparticles(0.0, 1e-32, sim.structure.DOS[mat], sim.structure.egrid)
            end
        else
            no_part = fill(get_thermalparticles(0.0,1e-32, sim.structure.DOS,sim.structure.egrid), sim.structure.dimension.length)
        end
        return no_part
    end
    if sim.electronictemperature.Enabled == true && sim.electronictemperature.AthermalElectron_ElectronCoupling == false
        if typeof(sim.structure.DOS) == Vector{spl} && sim.structure.Elemental_System == 1
            no_part = zeros(sim.structure.dimension.length)
            for j in eachindex(sim.structure.dimension.grid)
                no_part[j] = get_thermalparticles(0.0, 1e-32, sim.structure.DOS[j], sim.structure.egrid)
            end
        elseif typeof(sim.structure.DOS) == Vector{spl} && sim.structure.Elemental_System > 1
            no_part = zeros(sim.structure.dimension.length)
            for j in eachindex(sim.structure.dimension.grid)
                mat = mat_picker(sim.structure.dimension.grid[j], sim.structure.dimension.InterfaceHeight)
                no_part[j] = get_thermalparticles(0.0,1e-32, sim.structure.DOS[mat], sim.structure.egrid)
            end
        else
            no_part = fill(get_thermalparticles(0.0, 1e-32, sim.structure.DOS,sim.structure.egrid), sim.structure.dimension.length)
        end
        return no_part
    end
end

function write_electronictemperature(f,results,dim,mp,μs,sim,FD,relax)
    Tel = results["Tel"]
    Tph = results["Tph"]
    f["Temperature"] = Tel
    f["Distribution"] = FD
    if sim.Systems.ElectronTemperature == true
        f["dTel / dt"] = pp_derivative_Temp(Tel)
        f["Thermal Conductivity"] = pp_spatial_Tel(Tel,dim,Tph,mp)
        if sim.Interactions.ElectronPhonon == true
            f["Electron-Phonon Coupling"] = pp_electronphonon(Tel,Tph,cons,mp,μs,sim)
        end
        if sim.Interactions.ElectronElectron == false
            f["Electronic Heat Capacity"] = pp_electronheatcapacity(Tel,mp,cons,μs,sim)
            HDF5.API.h5l_create_soft("/Miscellaneous/Laser/Temporal Profile",fid_id,"/Electronic Temperature/Laser Source",HDF5.H5P_DEFAULT,HDF5.H5P_DEFAULT)
        elseif sim.Interactions.ElectronElectron == true
            f["Non-Eq Electron-Electron Energy Change"] =pp_neqelectronelectronenergychange(relax.ee,mp)
        end
    end
end

function pp_FermiDistribution(Tel,sim,μ)
    FD=zeros(length(Tel[:,1]),length(Tel[1,:]),length(sim.structure.egrid))
    Threads.@threads for i in eachindex(Tel[:,1])
        for j in eachindex(Tel[1,:])
            FD[i,j,:] .= FermiDirac(Tel[i,j],μ[i,j],sim.structure.egrid)
        end
    end
    return FD
end 

function pp_derivative_Temp(Temp)
    dTemp = zeros(size(Temp))
    for i in eachindex(dTemp[:,1])
        if i ==1
            dTemp[i,:] .= zeros(length(dTemp[i,:]))
        else
            dTemp[i,:] .= Temp[i,:] .- Temp[i-1,:]
        end
    end
    return dTemp
end

function pp_spatial_Tel(Tel,dim,Tph,mp)
    spat = zeros(size(Tel))
    Threads.@threads for i in eachindex(Tel[:,1])
        temp = zeros(length(Tel[1,:]))
        electrontemperature_conductivity(Tel[i,:],dim,Tph[i,:],mp,temp)
        spat[i,:] .= temp
    end
    return spat
end

function pp_electronphonon(Tel,Tph,cons,mp,μ,sim)
    eps = zeros(size(Tel))
    if sim.ParameterApprox.ElectronPhononCoupling == true
        Threads.@threads for i in eachindex(Tel[:,1])
            eps[i,:] .= nonlinear_electronphononcoupling.(cons.hbar,cons.kB,mp.λ,Ref(mp.DOS),Tel[i,:],μ[i,:],Tph[i,:])
        end
    else
        Threads.@threads for i in eachindex(Tel[:,1])
            eps[i,:] .= mp.g * (Tel[i,:].-Tph[i,:])
        end
    end
    return eps
end

function pp_electronheatcapacity(Tel,mp,cons,μ,sim)
    hc = zeros(size(Tel))
    if sim.ParameterApprox.ElectronHeatCapacity == true
        Threads.@threads for i in eachindex(Tel[:,1])
            hc[i,:] .= nonlinear_electronheatcapacity.(cons.kB,Tel[i,:],μ[i,:],Ref(mp.DOS))
        end
    else
        Threads.@threads for i in eachindex(Tel[:,1])
            hc[i,:] .= mp.γ *Tel[i,:]
        end
    end
    return hc
end

function pp_neqelectronelectronenergychange(relax,mp)
    uee = zeros(length(relax[1,1,:]),length(relax[1,:,1]))
    Threads.@threads for i in eachindex(uee[:,1])
        for j in eachindex(uee[i,:])
            uee[i,j] = elec_energychange(sim.structure.egrid,-1*relax[:,j,i],mp.DOS)
        end
    end
    return uee
end

function pp_chemicalpotential(Tel,n,sim)
    cp = zeros(size(Tel))
    if sim.athermalelectrons.Enabled == true
        Threads.@threads for i in eachindex(Tel[1,:])
            if typeof(sim.structure.DOS) == Vector{spl}
                cp[:,i] .= find_chemicalpotential.(n[:,i],Tel[:,i],Ref(sim.structure.DOS[i]),Ref(sim.structure.egrid))
            else
                cp[:,i] .= find_chemicalpotential.(n[:,i],Tel[:,i],Ref(sim.structure.DOS),Ref(sim.structure.egrid))
            end
        end
    else
        Threads.@threads for i in eachindex(Tel[1,:])
            if typeof(sim.structure.DOS) == Vector{spl} && sim.structure.Elemental_System == 1
                n = get_thermalparticles(0.0,1e-32,sim.structure.DOS[i],sim.structure.egrid)
                cp[:,i] .= find_chemicalpotential.(n,Tel[:,i],Ref(sim.structure.DOS[i]),Ref(sim.structure.egrid))
            elseif sim.structure.Elemental_System > 1
                idx = mat_picker(sim.structure.dimension.grid[i],sim.structure.dimension.InterfaceHeight)
                n = get_thermalparticles(0.0,1e-32,sim.structure.DOS[idx],sim.structure.egrid)
                cp[:,i] .= find_chemicalpotential.(n,Tel[:,i],Ref(sim.structure.DOS[idx]),Ref(sim.structure.egrid))
            else
                n = get_thermalparticles(0.0,1e-32,sim.structure.DOS,sim.structure.egrid)
                cp[:,i] .= find_chemicalpotential.(n,Tel[:,i],Ref(sim.structure.DOS),Ref(sim.structure.egrid))
            end
        end
    end
    return cp
end

function write_phononictemperature(f,results,sim,mp,fid_id)
    Tph = results["Tph"]
    f["Temperature"] = Tph
    if sim.Systems.PhononTemperature == true
        f["dTph / dt"] = pp_derivative_Temp(Tph)
        HDF5.API.h5l_create_soft("/Electronic Temperature/Electron-Phonon Coupling",fid_id,"/Phononic Temperature/Electron-Phonon Coupling",HDF5.H5P_DEFAULT,HDF5.H5P_DEFAULT)
        f["Phononic Heat Capacity"] = pp_phononicheatcapacity(Tph,mp,results["n"],cons)
        if sim.Systems.NonEqElectrons == true
            f["Non-Eq Electron-Phonon Energy Change"] = pp_neqelectronphononenergychange(results["fneq"],mp)
        end
    end
end

function pp_phononicheatcapacity(Tph,mp,n,cons)
    if sim.ParameterApprox.PhononHeatCapacity == true
        hc = zeros(size(Tph))
        if sim.Systems.NonEqElectrons == true
            Threads.@threads for i in eachindex(Tph[:,1])
                hc[i,:] .= nonlinear_phononheatcapacity.(Tph[i,:],n[i,:],cons.kB,mp.θ)
            end
        else
            Threads.@threads for i in eachindex(Tel[:,1])
                hc[i,:] .= nonlinear_phononheatcapacity.(Tph[i,:],n,cons.kB,mp.θ)
            end
        end
        return hc
    else
        return mp.Cph
    end
end

function pp_neqelectronphononenergychange(fneq,mp)
    uep = zeros(length(fneq[1,1,:]),length(fneq[1,:,1])) 
    Threads.@threads for i in eachindex(fneq[1,1,:])
        for j in eachindex(fneq[1,:,1])
            uep[i,j] = neq_electronphonontransfer(fneq[:,j,i],sim.structure.egrid,mp.τep,mp.DOS)
        end
    end
    return uep
end

function write_numberofparticles(f,results,sim,mp,relax,μ)
    if sim.Systems.NonEqElectrons == true
        f["Number of Electrons"] = results["noe"]
        if sim.Interactions.ElectronElectron == true
            f["dn / dt"] = pp_changeinparticles(relax.ee,mp,μ)
        end
    else
        f["Number of Electrons"] = mp.n0
    end
end

function pp_changeinparticles(relax,mp,μ)
    Δn = zeros(size(μ))
    Threads.@threads for i in eachindex(relax[1,1,:])
        temp = zeros(size(Δn[i,:]))
        noe_func(-1*relax[:,:,i],μ[i,:],mp,temp)
        Δn[i,:] .= temp
    end
    return Δn
end

function write_noneqelectrons(f,results,sim,FD,μ,mp,cons,las)
    fneq = results["fneq"]
    ftot = fneq .+ FD
    f["Distribution"] = fneq
    f["Total Distribution"] = ftot
    f["dfneq / dt"] = pp_dfneqdt(fneq)
    f["Excitation"] = pp_athemexcitation(ftot,mp,las)
    if sim.Interactions.ElectronElectron == true
        elel =  pp_athemelectronelectronrelaxation(results["Tel"],fneq,results["n"],μ,mp,cons)
        f["Electron-Electron Relaxation"] = elel
        relax = (;relax...,ee=elel)
    end
    if sim.Interactions.ElectronPhonon == true
        elph = pp_athemelectronphononrelaxation(fneq,mp)
        f["Electron-Phonon Relaxation"] = elph
        relax=(;relax...,ep=elph)
    end
    return relax
end

function pp_dfneqdt(fneq)
    dfneq = zeros(size(fneq))
    Threads.@threads for i in eachindex(fneq[1,1,:])
        if i == 1
            dfneq[:,:,i] .= dfneq[:,:,i]
        else
            @. dfneq = fneq[:,:,i] - fneq[:,:,i-1]
        end
    end
    return dfneq
end

function pp_athemexcitation(ftot,mp,las)
    excite = zeros(size(ftot))
    Threads.@threads for i in eachindex(ftot[1,1,:])
        for j in eachindex(ftot[1,:,1])
            excite[:,j,i] .= athemexcitation(ftot[:,j,i],sim.structure.egrid,mp.DOS,las.hv)
        end
    end
    return FD
end 

function pp_athemelectronelectronrelaxation(Tel,fneq,n,μ,mp,cons)
    rel = zeros(size(fneq))
    Threads.@threads for i in eachindex(fneq[1,1,:])
        temp = zeros(size(rel[:,:,i]))
        relax_func(Tel[i,:],fneq[:,:,i],n[i,:],μ[i,:],mp,cons,temp)
        rel[:,:,i] .= -1*temp
    end
    return rel
end

function pp_athemelectronphononrelaxation(fneq,mp)
    rel = zeros(size(fneq))
    Threads.@threads for i in eachindex(fneq[1,1,:])
        rel[:,:,i] .= -fneq[:,:,i]./mp.τep
    end
    return rel
end

function pp_temporalprofile(las,dim,timepoints,mp)
    temp_prof = zeros(length(timepoints),length(dim.grid))
    if typeof(las) == Gaussian
        for i in eachindex(timepoints)
            temp_prof[i,:] = sqrt(4*log(2)/pi)/las.FWHM*exp(-4*log(2)*timepoints[i]^2/las.FWHM^2)*(1-mp.R)*las.Power*1/(mp.ϵ*(1-exp(-dim.grid[end]/mp.ϵ))).*exp.(-dim.grid./mp.ϵ)
        end
    end
    return temp_prof
end

function write_minimum(f,results,FD,sim)
    write_dataset(f["Electronic Temperature"],"Temperature",results["Tel"])
    #f["Electronic Temperature"]["Temperature"][:] = results["Tel"]
    write_dataset(f["Electronic Temperature"],"Distribution",FD)
    #f["Electronic Temperature"]["Distribution"][:] = FD
    write_dataset(f["Phononic Temperature"],"Temperature",results["Tph"])
    #f["Phononic Temperature"]["Temperature"][:] = results["Tph"]
    write_dataset(f,"Particle Number",results["noe"])
    #f["Particle Number"][:] = results["noe"]
    write_dataset(f["Athermal Electrons"],"Non-Equilibrium Distribution",results["fneq"])
    #f["Athermal Electrons"]["Non-Equilibrium Distribution"][:] = results["fneq"]
    if sim.athermalelectrons.Enabled == true
        if size(FD,1) != size(results["fneq"],1)
            FD = repeat(FD,size(results["fneq"],1), 1, 1)
            FD_tot = results["fneq"].+FD
            write_dataset(f["Electronic Distribution"],"Total Distribution",FD_tot)
            #f["Electronic Distribution"]["Total Distribution"][:] = FD_tot
        else
            FD_tot = results["fneq"].+FD
            write_dataset(f["Electronic Distribution"],"Total Distribution",FD_tot)
            #f["Electronic Distribution"]["Total Distribution"][:] = FD_tot
        end
    end
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

    
function output_functions(f,sim)
    sys = function_builder(sim)
    for i in keys(sys)
        f[i] = string(sys[i])
    end
    simulation_expr = simulation_construction(sys,sim)
    f["Total Equation Block"] = string(simulation_expr)
end
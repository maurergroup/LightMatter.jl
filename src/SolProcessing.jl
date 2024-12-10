function post_production(sol,file_name,initial_temps,output,sim,mp,las,dim,cons)
    fid = create_datafile_and_structure(file_name)
    fid_id = fid.id
    write_simsettings(fid["Settings"],sim)
    write_materialproperties(fid["Parameters"],mp,dim,cons)

    results = seperate_results(sol,initial_temps,mp,sim)
    fid["Miscellaneous"]["Time Span"] = results["times"]
    μs = pp_chemicalpotential(results["Tel"],results["n"],mp,cons)
    fid["Miscellaneous"]["Chemical Potential"] = μs
    FD = pp_FermiDistribution(results["Tel"],mp,cons,μs)
    write_laser(fid["Miscellaneous"],las,dim,results["times"],mp)

    create_group(fid["Miscellaneous"],"System Equations")
    output_functions(fid["Miscellaneous"]["System Equations"],sim,las,dim)

    if output=="full"
        relax=()
        if sim.Systems.NonEqElectrons == true
            relax=write_noneqelectrons(fid["Non Eq Electrons"],results,sim,FD,μs,mp,cons,las)
        end
        write_electronictemperature(fid["Electronic Temperature"],results,dim,mp,μs,sim,FD,relax)
        write_phononictemperature(fid["Phononic Temperature"],results,sim,mp,fid_id)
        write_numberofparticles(fid["Number Of Particles"],results,sim,mp,relax,μs)
        close(fid)
    elseif output=="minimum"
        write_minimum(fid,results,FD,sim)
        close(fid)
    else 
        println("The only output options are full and minimum. Files cannot be overwritten so the
        file has been deleted")
        close(fid)
        rm(file_name)
    end
end

function create_datafile_and_structure(file_name)
    fid = h5open(file_name,"w")
    create_group(fid,"Settings")
    create_group(fid,"Parameters")
    create_group(fid,"Electronic Temperature")
    create_group(fid,"Phononic Temperature")
    create_group(fid,"Number Of Particles")
    create_group(fid,"Non Eq Electrons")
    create_group(fid,"Miscellaneous")
    return fid
end

function dict_to_hdf5(g,d)
    for (key,value) in d
        g[key] = value
    end
end

function write_simsettings(f,sim)
    Ints = Dict{String,Union{Float64,Bool}}(String(key)=>convert(Bool,getfield(sim.Interactions, key)) for key ∈ fieldnames(Interaction))
    create_group(f,"Interactions")
    dict_to_hdf5(f["Interactions"],Ints)
    ParamApprox = Dict{String,Union{Float64,Bool}}(String(key)=>convert(Bool,getfield(sim.ParameterApprox, key)) for key ∈ fieldnames(ParameterApproximation))
    create_group(f,"Parameter Approximations")
    dict_to_hdf5(f["Parameter Approximations"],ParamApprox)
    comps = Dict{String,Union{Float64,Bool}}(String(key)=>convert(Bool,getfield(sim.Systems, key)) for key ∈ fieldnames(SystemComponents))
    create_group(f,"System Components")
    dict_to_hdf5(f["System Components"],comps)
end

function write_materialproperties(f,mp,dim,cons)
    tuple_mp = [("Extinction Coefficient / nm^-1",mp.ϵ),
                ("Valence Band Minimum / eV",mp.FE),
                ("Electronic Specific Heat Capacity / eV/fs/nm^3/K^2",mp.γ),
                ("Debye Temperature / K",mp.θ),
                ("Number of Atoms / nm^3",mp.n),
                ("Room Temperature Diffusive Thermal Conductivity / eV/fs/m/K",mp.κ),
                ("Number of Free Electrons per atom",mp.ne),
                ("Effective Mass of Electrons",mp.effmass),
                ("Second Moment of Phonon Spectrum × Electron-Phonon Mass Enhancement Factor / eV^2",mp.λ),
                ("Linear Electron-Phonon Coupling Constant / eV/fs/nm^3/K",mp.g),
                ("Material Dependent Fermi Liquid Lifetime scalar / fs^-1",mp.τ),
                ("Electron-Phonon Relaxation Time / fs^-1",mp.τep),
                ("Energy Grid for solving Distributions",mp.egrid),
                ("Constant Phononic Heat Capacity",mp.Cph)]
    matpat = Dict(tuple_mp)
    create_group(f,"Material Parameters")
    dict_to_hdf5(f["Material Parameters"],matpat)
    constants = Dict("Boltzmann Constant / eV/K"=>cons.kB,"Reduced Planck's Constant / eVfs"=>cons.hbar)
    create_group(f,"Constants")
    dict_to_hdf5(f["Constants"],constants)
    dimensions = write_dimsdict(dim)
    create_group(f,"Dimensions")
    dict_to_hdf5(f["Dimensions"],dimensions)
    write_DOS(f,mp)
end

function write_DOS(f,mp)
    create_group(f,"Density Of States")
    egrid = collect(range(-10,10,step=0.01))
    DOS = zeros(length(mp.DOS),length(egrid))
    for i in eachindex(mp.DOS)
        DOS[i,:] = mp.DOS[i](egrid)
    end
    f["Density Of States"]["DOS"] = DOS
    f["Density Of States"]["Energy Axis"] = egrid

end

function write_dimsdict(dim)
    if typeof(dim) == Homogeneous
        return Dict("Dimension" => 0)
    elseif typeof(dim) == Linear
        return Dict("Dimension" => 1, "Number of Points" => dim.length, "Spacing" => dim.dz, "Grid" => dim.grid)
    end
end

function seperate_results(sol,initial_temps,mp,sim)
    l = length(sol[:,1].x)
    vals =  Dict{String,Any}("Tel"=>0.0,"Tph"=>0.0,"fneq"=>0.0,"n"=>0.0,"times"=>sol.t)
    if l == 1
        vals["Tel"] = [initial_temps["Tel"]]
        vals["Tph"] = [initial_temps["Tel"]]
        vals["fneq"] = cat(getindex.(getfield.(sol.u, :x), 1)...,dims=3)
        vals["n"] = mp.n0
    elseif l== 2
        vals["Tel"] = stack(getindex.(getfield.(sol.u, :x), 1),dims=1)
        vals["Tph"] = stack(getindex.(getfield.(sol.u, :x), 2),dims=1)
        vals["n"] = mp.n0
    elseif l==3
        vals["Tel"] = stack(getindex.(getfield.(sol.u, :x), 2),dims=1)
        vals["Tph"] = [initial_temps["Tel"]]
        vals["fneq"] = cat(getindex.(getfield.(sol.u, :x), 3)...,dims=3)
        vals["n"] = stack(getindex.(getfield.(sol.u, :x), 1),dims=1)
    elseif l==4
        if sim.ParameterApprox.EmbeddingMethod == true
            vals["Tel"] = stack(getindex.(getfield.(sol.u, :x), 3),dims=1)
            vals["Tph"] = stack(getindex.(getfield.(sol.u, :x), 4),dims=1)
            vals["fneq"] = cat(getindex.(getfield.(sol.u, :x), 2)...,dims=3)
            vals["n"] = stack(getindex.(getfield.(sol.u, :x), 1),dims=1)
        else
            vals["Tel"] = stack(getindex.(getfield.(sol.u, :x), 1),dims=1)
            vals["Tph"] = stack(getindex.(getfield.(sol.u, :x), 3),dims=1)
            vals["fneq"] = cat(getindex.(getfield.(sol.u, :x), 4)...,dims=3)
            vals["n"] = stack(getindex.(getfield.(sol.u, :x), 2),dims=1)
        end
    end
    return vals
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

function pp_FermiDistribution(Tel,mp,cons,μ)
    FD=zeros(length(Tel[1,:]),length(mp.egrid),length(Tel[:,1]))
    Threads.@threads for i in eachindex(Tel[:,1])
        for j in eachindex(Tel[1,:])
            FD[j,:,i] .= FermiDirac(Tel[i,j],μ[i,j],cons.kB,mp.egrid)
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
            uee[i,j] = elec_energychange(mp.egrid,-1*relax[:,j,i],mp.DOS,mp.FE)
        end
    end
    return uee
end

function pp_chemicalpotential(Tel,n,mp,cons)
    cp = zeros(size(Tel))
    if typeof(n) == Float64
        Threads.@threads for i in eachindex(Tel[1,:])
            cp[:,i] .= find_chemicalpotential.(n,Tel[:,i],Ref(mp.DOS[i]),cons.kB,mp.FE)
        end
    else
        Threads.@threads for i in eachindex(Tel[1,:])
            cp[:,i] .= find_chemicalpotential.(n[:,i],Tel[:,i],Ref(mp.DOS[i]),cons.kB,mp.FE)
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
            uep[i,j] = neq_electronphonontransfer(fneq[:,j,i],mp.egrid,mp.τep,mp.FE,mp.DOS)
        end
    end
    return uep
end

function write_numberofparticles(f,results,sim,mp,relax,μ)
    if sim.Systems.NonEqElectrons == true
        f["Number of Electrons"] = results["n"]
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
            excite[:,j,i] .= athemexcitation(ftot[:,j,i],mp.egrid,mp.DOS,las.hv,mp.FE)
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

function write_laser(f,las,dim,timepoints,mp)
    create_group(f,"Laser")
    tuple_las = [("Laser Form",String(Symbol(typeof(las)))),
                 ("Full-Width at Half Maximum / fs",las.FWHM),
                 ("Fluence / eV/nm^2",las.Power),
                 ("Photon Energy / eV",las.hv),
                 ("Reflectivity",las.R),
                 ("Transport Type",las.Transport)]
    las_dict=Dict(tuple_las)
    dict_to_hdf5(f["Laser"],las_dict)
    laser=laser_factory(las,dim)
    f["Laser"]["Temporal Profile"] = pp_temporalprofile(las,dim,timepoints,mp)
end

function pp_temporalprofile(las,dim,timepoints,mp)
    temp_prof = zeros(length(timepoints),length(dim.grid))
    if typeof(las) == Gaussian
        for i in eachindex(timepoints)
            temp_prof[i,:] = sqrt(4*log(2)/pi)/las.FWHM*exp(-4*log(2)*timepoints[i]^2/las.FWHM^2)*(1-las.R)*las.Power*1/(mp.ϵ*(1-exp(-dim.grid[end]/mp.ϵ))).*exp.(-dim.grid./mp.ϵ)
        end
    end
    return temp_prof
end

function write_minimum(f,results,FD,sim)
    f["Electronic Temperature"]["Temperature"] = results["Tel"]
    f["Electronic Temperature"]["Distribution"] = FD
    f["Phononic Temperature"]["Temperature"] = results["Tph"]
    f["Number Of Particles"]["Particle Number"] = results["n"]
    f["Non Eq Electrons"]["Non-Equilibrium Distribution"] = results["fneq"]
    if sim.Systems.NonEqElectrons == true
        f["Non Eq Electrons"]["Total Distribution"] = results["fneq"].+FD
    end
end
    
function output_functions(f,sim,las,dim)
    sys,key_list = function_builder(sim,las,dim,sys_output=true)
    f["Name of Method"] = method_name(key_list)
    for i in key_list
        f[i] = string(sys)
    end
end

function method_name(key_list)
    if length(key_list) == 1
        return "Athermal Electron Generation"
    elseif length(key_list) == 2
        return "Two-Temperature Model"
    elseif length(key_list) == 3
        return "Athermal Electrons with Electron Relaxation"
    elseif length(key_list) == 4
        return "AthEM with Electron & Phonon Relaxation"
    elseif length(key_list) == 7
        return "Embedded AthEM inside Two-Temperature Model"
    end
end
function electrontemperature_factory(sim::SimulationSettings,laser::Expr,dim::Dimension)
    if sim.Systems.NonEqElectrons == false
        HeatCapacity = electrontemperature_heatcapacity(sim)
        ElecPhon = electronphonon_coupling(sim)
        Spatial = :(cond) 
        return build_electronTTM(laser,Spatial,ElecPhon,HeatCapacity)
    elseif sim.Systems.NonEqElectrons == true
        Δu = athem_electempenergychange(sim,dim)
        return build_athemelectron(Δu) 
    end
end

function build_electronTTM(Source,Spatial,ElecPhon,HeatCapacity)
    return Expr(:call,:/,Expr(:call,:+,Source,Spatial,ElecPhon),HeatCapacity)
end

function electrontemperature_heatcapacity(sim::SimulationSettings)
    if sim.ParameterApprox.ElectronHeatCapacity == true
        return :(nonlinear_electronheatcapacity(cons.kB,Tel,μ,DOS,mp.egrid))
    else
        return :(mp.γ*Tel)
    end
end

function nonlinear_electronheatcapacity(kB,Tel,μ,DOS,egrid)
    return extended_Bode(dFDdT(kB,Tel,μ,egrid).*DOS(egrid).*egrid,egrid)
end

function electronphonon_coupling(sim)
    if sim.Interactions.ElectronPhonon == true
        if sim.ParameterApprox.ElectronPhononCoupling==true
            return :(nonlinear_electronphononcoupling(cons.hbar,cons.kB,mp.λ,DOS,Tel,μ,Tph,mp.egrid))
        else
            return :(-mp.g*(Tel-Tph))
        end
    else
        return 0.0
    end
end

function nonlinear_electronphononcoupling(hbar,kB,λ,DOS,Tel,μ,Tph,egrid)
    prefac=pi*kB*λ/DOS(μ)/hbar
    g=prefac.*extended_Bode(DOS(egrid).^2 .*-dFDdE(kB,Tel,μ,egrid),egrid)
    return -g*(Tel-Tph)
end

function build_athemelectron(Δu)
    return :( 1/(c_T(μ,Tel,DOS,cons.kB,mp.egrid)*p_μ(μ,Tel,DOS,cons.kB,mp.egrid)-p_T(μ,Tel,DOS,cons.kB,mp.egrid)*c_μ(μ,Tel,DOS,cons.kB,mp.egrid))*(p_μ(μ,Tel,DOS,cons.kB,mp.egrid)*$Δu-c_μ(μ,Tel,DOS,cons.kB,mp.egrid)*Δn))
end

function elec_energychange(egrid,relax_dis,DOS)
    return get_internalenergy(relax_dis,DOS,egrid)
end

function athem_electempenergychange(sim,dim)
    args = Vector{Union{Expr,Symbol}}(undef,0)
    push!(args,:(elec_energychange(mp.egrid,relax_dis,DOS)))
    if sim.Systems.PhononTemperature == true
       push!(args,:(nonlinear_electronphononcoupling(cons.hbar,cons.kB,mp.λ,DOS,Tel,μ,Tph,mp.egrid)))
    end
    if typeof(dim) != Homogeneous
        push!(args,:(cond))
    end
    return :(.+($(args...)))
end

function electrontemperature_conductivity(Tel,dim::Dimension,Tph,mp,cond)
    Depthderivative(Tel,dim.dz,cond)
    cond[1]=0.0
    cond[end]=0.0
    K=mp.κ*Tel./Tph
    Depthderivative((cond.*K),dim.dz,cond)
end

function electrontemperature_conductivity(Tel,dim::Homogeneous,Tph,mp,cond)
    return [0.0]
end

function Depthderivative(vec,dz,Diff)
    for i in 2:length(vec)-1
        Diff[i]=(vec[i+1]-vec[i-1])/(2*dz)
    end
    Diff[1] = (vec[2]-vec[1])/dz
    Diff[end] = (vec[end]-vec[end-1])/dz
end

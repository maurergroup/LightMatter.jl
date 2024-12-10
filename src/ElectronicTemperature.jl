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
        return :(nonlinear_electronheatcapacity(cons.kB,Tel,μ,DOS))
    else
        return :(mp.γ*Tel)
    end
end

function nonlinear_electronheatcapacity(kB,Tel,μ,DOS)
    int(u,p) = dFDdT(kB,Tel,μ,u)*DOS(u)*u
    prob = IntegralProblem(int,(μ-(60*Tel/10000),μ+(60*Tel/10000)))
    return solve(prob,HCubatureJL(initdiv=50);reltol=1e-5,abstol=1e-5).u
end

function electronphonon_coupling(sim)
    if sim.Interactions.ElectronPhonon == true
        if sim.ParameterApprox.ElectronPhononCoupling==true
            return :(nonlinear_electronphononcoupling(cons.hbar,cons.kB,mp.λ,DOS,Tel,μ,Tph))
        else
            return :(-mp.g*(Tel-Tph))
        end
    else
        return 0.0
    end
end

function nonlinear_electronphononcoupling(hbar,kB,λ,DOS,Tel,μ,Tph)
    prefac=pi*kB*λ/DOS(μ)/hbar
    int(u,p) = DOS(u)^2*-dFDdE(kB,Tel,μ,u)
    prob = IntegralProblem(int,(μ-(10*Tel/10000),μ+(10*Tel/10000)))
    g=prefac.*solve(prob,HCubatureJL(initdiv=50);reltol=1e-5,abstol=1e-5).u
    return -g*(Tel-Tph)
end

function build_athemelectron(Δu)
    return :( 1/(c_T(μ,Tel,DOS,cons.kB)*p_μ(μ,Tel,DOS,cons.kB)-p_T(μ,Tel,DOS,cons.kB)*c_μ(μ,Tel,DOS,cons.kB))*(p_μ(μ,Tel,DOS,cons.kB)*$Δu-c_μ(μ,Tel,DOS,cons.kB)*Δn))
end

function elec_energychange(egrid,relax_dis,DOS,FE)
    spl = get_interpolate(egrid,relax_dis)
    return get_internalenergyspl(spl,DOS,FE)
end

function athem_electempenergychange(sim,dim)
    args = Vector{Union{Expr,Symbol}}(undef,0)
    push!(args,:(elec_energychange(mp.egrid,relax_dis,DOS,mp.FE)))
    if sim.Systems.PhononTemperature == true
       push!(args,:(nonlinear_electronphononcoupling(cons.hbar,cons.kB,mp.λ,DOS,Tel,μ,Tph)))
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

function embedded_AthEM_conductivity(Tel,Tph,dim,mp)
    dTel = (Tel[3].-Tel[1])./dim.dz
    K=mp.κ*Tel[2]./Tph[2]
    dQ = dTel .* K
    return (dQ-0.0)/dim.dz
end


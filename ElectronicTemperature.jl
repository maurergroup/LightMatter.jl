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
        return :(nonlinear_electronheatcapacity(cons.kB,Tel,μ,mp.DOS))
    else
        return :(mp.γ*Tel)
    end
end

function nonlinear_electronheatcapacity(kB::Real,Tel::Real,μ::Real,DOS::spl)
    int(u,p) = dFDdT(kB,Tel,μ,u)*DOS(u)*u
    prob = IntegralProblem(int,(μ-(60*Tel/10000),μ+(60*Tel/10000)))
    return solve(prob,HCubatureJL(initdiv=2);reltol=1e-5,abstol=1e-5).u
end

function electronphonon_coupling(sim)
    if sim.Interactions.ElectronPhonon == true
        if sim.ParameterApprox.ElectronPhononCoupling==true
            return :(nonlinear_electronphononcoupling(cons.hbar,cons.kB,mp.λ,mp.DOS,Tel,μ,Tph))
        else
            return :(-mp.g*(Tel-Tph))
        end
    else
        return 0.0
    end
end

function nonlinear_electronphononcoupling(hbar::Real,kB::Real,λ::Real,DOS::spl,Tel::Real,μ::Real,Tph::Real)
    prefac=pi*kB*λ/DOS(μ)/hbar
    int(u,p) = DOS(u)^2*-dFDdE(kB,Tel,μ,u)
    prob = IntegralProblem(int,(μ-(60*Tel/10000),μ+(60*Tel/10000)))
    g=prefac.*solve(prob,HCubatureJL(initdiv=2);reltol=1e-5,abstol=1e-5).u
    return -g*(Tel-Tph)
end

function electrontemperature_conductivity(dim::Dimension)
    if typeof(dim) == Homogeneous
        return 0.0
    else 
        return :(electronicthermalconductivity())
    end
end

function build_athemelectron(Δu)
    return :( 1/(c_T(μ,Tel,mp.DOS,cons.kB)*p_μ(μ,Tel,mp.DOS,cons.kB)-p_T(μ,Tel,mp.DOS,cons.kB)*c_μ(μ,Tel,mp.DOS,cons.kB))*(p_μ(μ,Tel,mp.DOS,cons.kB)*$Δu-c_μ(μ,Tel,mp.DOS,cons.kB)*Δn))
end

function elec_energychange(egrid,relax_dis,μ,DOS,u0,FE)
    spl = get_interpolate(egrid,relax_dis)
    return get_internalenergyspl(μ,spl,DOS,u0,FE)
end

function athem_electempenergychange(sim,dim)
    args = Vector{Union{Expr,Symbol}}(undef,0)
    push!(args,:(elec_energychange(mp.egrid,relax_dis,μ,mp.DOS,mp.u0,mp.FE)))
    if sim.Systems.PhononTemperature == true
       push!(args,:(nonlinear_electronphononcoupling(cons.hbar,cons.kB,mp.λ,mp.DOS,Tel,μ,Tph)))
    end
    if typeof(dim) != Homogeneous
        push!(args,:(cond))
    end
    return :(+($(args...)))
end

function electrontemperature_conductivity(Tel,dim,Tph,mp)
    Tel_spl = get_interpolate(dim.grid,Tel)
    dTeldz = DataInterpolations.derivative.(Ref(Tel_spl),dim.grid)
    K=mp.κ*(Tel./Tph)
    Q_spl=get_interpolate(dim.grid,dTeldz.*K)
    cond = DataInterpolations.derivative.(Ref(Q_spl),dim.grid)
    cond[end]=0.0
    return cond
end

function electrontemperature_conductivity(Tel,dim::Homogeneous,Tph,mp)
    return [0.0]
end
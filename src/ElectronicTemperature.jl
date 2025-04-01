"""
    electrontemperature_factory(sim::Simulation,laser::Expr)
    This function takes the Simulation struct and a laser expression and returns an assembled
    expression for how the electronic temperature should evolve over time. 
"""
function electrontemperature_factory(sim::Simulation,laser::Expr)
    if sim.electronictemperature.AthermalElectron_ElectronCoupling == false
        HeatCapacity = electrontemperature_heatcapacity(sim)
        ElecPhon = electronphonon_coupling(sim)
        return build_electronTTM(sim,laser,ElecPhon,HeatCapacity)
    elseif sim.electronictemperature.AthermalElectron_ElectronCoupling == true
        Δu = athem_electempenergychange(sim)
        return build_athemelectron(Δu) 
    end
end
"""
    build_electronTTM(Source::Expr,Spatial::Symbol,ElecPhon::Expr,HeatCapacity::Expr)
    Builds a combine expression of the expressions for the energy input, electorn-phonon coupling,
    thermal conductivity and the electronic heat capacity. The first three terms are summed to find
    the change in internal energy of the thermal system and divided by the heat capacity to change the 
    internal energy into a temperature.
"""
function build_electronTTM(sim::Simulation,Source::Expr,ElecPhon::Expr,HeatCapacity::Expr)
    args = Union{Expr,Symbol,Real}[Source,ElecPhon]
    if sim.electronictemperature.Conductivity == true
        push!(args,:Tel_cond)
    end
    return Expr(:call,:/,Expr(:call,:+,args...),HeatCapacity)
end
"""
    electrontemperature_heatcapacity(sim::Simulation)
    Returns the expression for how the electronic temperatures heat capacity should be 
    calculated. This can be done either via a linear relation between the speciic heat(γ) and
    Tel or via a more accurate but complex non-linear relationship. The keyword to be set 
    in define_simulation_settings is nlelecheat = true.
"""
function electrontemperature_heatcapacity(sim::Simulation)
    if sim.electronictemperature.ElectronicHeatCapacity == :nonlinear
        return :(Lightmatter.nonlinear_electronheatcapacity(Tel,μ,DOS,sim.structure.egrid))
    elseif sim.electronictemperature.ElectronicHeatCapacity == :linear
        return :(sim.electronictemperature.γ*Tel)
    end
end
"""
    nonlinear_electronheatcapacity(kB::Real,Tel::Real,μ::Real,DOS::spl,egrid::Vector{<:Real})
    This is the function that returns the nonlinear variant of the electronic heat capacity.
    It is defined as the internal energy of the derivative of the Fermi distribution with respect
    to temperature.
"""
function nonlinear_electronheatcapacity(Tel::Real,μ::Real,DOS::spl,egrid::Vector{<:Real})
    return extended_Bode(dFDdT(Tel,μ,egrid).*DOS(egrid).*egrid,egrid)
end
"""
    electronphonon_coupling(sim::Simulation)
    Returns the expression for how the electron-phonon coupling should be evaulated. It can use a constant
    electron-phonon parameter or a non-linear expression. The keyword to be set in define_simulation_settings
    is nlelecphon = true.
"""
function electronphonon_coupling(sim::Simulation)
    if sim.electronictemperature.Electron_PhononCoupling == true
        if sim.electronictemperature.ElectronPhononCouplingValue == :variable
            return :(Lightmatter.nonlinear_electronphononcoupling(sim.electronictemperature.λ,sim.electronictemperature.ω,DOS,Tel,μ,Tph,sim.structure.egrid))
        elseif sim.electronictemperature.ElectronPhononCouplingValue == :constant
            return :(-sim.electronictemperature.g*(Tel-Tph))
        end
    else
        return 0.0
    end
end
"""
    nonlinear_electronphononcoupling(λ::Real,ω::Real,DOS::spl,Tel::Real,μ::Real,Tph::Real,egrid::Vector{<:Real})
    Calculates the non-linear electron phonon coupling by calculating g. The equation for which can be found as equation 8 in
    Z. Lin, L. V. Zhigilei and V. Celli, Phys. Rev. B, 2008, 77, 075133.
"""
function nonlinear_electronphononcoupling(λ::Real,ω::Real,DOS::spl,Tel::Real,μ::Real,Tph::Real,egrid::Vector{<:Real})
    prefac=pi*Constants.kB*λ*ω/DOS(μ)/Constants.ħ
    g=prefac.*extended_Bode(DOS(egrid).^2 .*-dFDdE(Tel,μ,egrid),egrid)
    return -g*(Tel-Tph)
end
"""
    build_athemelectron(Δu::Expr)
    This function returns the expression for how the temperature changes during AthEM, this is different to the TTM because
    both the internal energy and number of particles changes so both must be accounted for. This leads to the rather complex 
    expression used here. 
"""
function build_athemelectron(Δu::Expr)
    return :( 1/(Lightmatter.c_T(μ,Tel,DOS,sim.structure.egrid)*Lightmatter.p_μ(μ,Tel,DOS,sim.structure.egrid)
    -Lightmatter.p_T(μ,Tel,DOS,sim.structure.egrid)*Lightmatter.c_μ(μ,Tel,DOS,sim.structure.egrid))*
    (Lightmatter.p_μ(μ,Tel,DOS,sim.structure.egrid)*$Δu-Lightmatter.c_μ(μ,Tel,DOS,sim.structure.egrid)*Δn))
end
"""
    elec_energychange(egrid::Vector{<:Real},relax_dis::Vector{<:Real},DOS::spl)
    Calculates the internal energy gained from the relaxation of the athermal electrons into the thermal bath by finding the
    internal energy of relax_dis.
"""
function elec_energychange(egrid::Vector{<:Real},relax_dis::Vector{<:Real},DOS::spl)
    return Lightmatter.get_internalenergy(relax_dis,DOS,egrid)
end
"""
    athem_electempenergychange(sim::Simulation)
    Creates the expression for how the internal energy of the electron temperature should increase during an AthEM calculation.
    This includes the source term from AthEM as well as thermal conductivity and electron-phonon coupling. 
"""
function athem_electempenergychange(sim::Simulation)
    args = Vector{Union{Expr,Symbol}}(undef,0)
    push!(args,:(Lightmatter.elec_energychange(sim.structure.egrid,relax_dis,DOS)))
    if sim.Systems.PhononTemperature == true
       push!(args,electronphonon_coupling(sim::Simulation))
    end
    if sim.electronictemperature.Conductivity == true
        push!(args,:(cond))
    end
    return Expr(:call,:+,args...)
end
"""
    electrontemperature_conductivity!(Tel::Vector{<:Real},κ::Real,dz::Real,Tph::Vector{<:Real},cond::Vector{<:Real})
    Calculates the thermal energy passing further into a bulk slab due to thermal conductivity of the electronic bath. 
    The derivative of the temperature with respect to distance is set to 0 at the boundaries.
"""
function electrontemperature_conductivity!(Tel::Vector{<:Real},κ::Real,dz::Real,Tph::Vector{<:Real},cond::Vector{<:Real})
    depthderivative!(Tel,dz,cond)
    cond[1] = 0.0
    cond[end] = 0.0
    K=κ.*Tel./Tph
    depthderivative!((cond.*K),dz,cond)
end
"""
    depthderivative!(vec::Vector{<:Real},dz::Real,Diff::Vector{<:Real})
    Calculates the finite differences derivative of the electronic temperature due to position. This uses central difference
    for all elements other than the first and last which use forward and reverse difference respectively to remove the need for
    environmental boundary conditions. 
"""
function depthderivative!(vec::Vector{<:Real},dz::Real,Diff::Vector{<:Real})
    for i in 2:length(vec)-1
        Diff[i]=(vec[i+1]-vec[i-1])/(2*dz)
    end
    Diff[1] = (vec[2]-vec[1])/dz
    Diff[end] = (vec[end]-vec[end-1])/dz
end
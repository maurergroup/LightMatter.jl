"""
    FermiDirac(Tel::Real, μ::Real, E::Union{Vector{<:Real},Real}) 
    
    Returns the thermal occupation of electrons at a given temperature, chemical potential and energy.
    If a vector of energies is given then it will return the distribution across that range.

    # Arguments
    - 'Tel': Thermal electronic temperature
    - 'μ': Chemical Potential
    - 'E': Energy value or range

    # Returns
    - Value or vector of the Fermi-Dirac distribution
"""
@inline FermiDirac(Tel::Real, μ::Real, E::Union{Vector{<:Real},Real}) = 1 ./ (exp.((E.-μ) ./ (Constants.kB*Tel)).+1)
"""
    dFDdE(Tel::Real, μ::Real, E::Union{Vector{<:Real},Real})
    
    Returns the change in the Fermi distribution with respect to energy at the given energy value or range.

    # Arguments
    - 'Tel': Thermal electronic temperature
    - 'μ': Chemical Potential
    - 'E': Energy value or range

    # Returns
    - Value or vector of the partial derivative of the Fermi distribution with respect to energy
"""
@inline function dFDdE(Tel::Real, μ::Real, E::Union{Vector{<:Real},Real})
    return -exp.((E.-μ)./(Constants.kB*Tel)) ./ (Constants.kB*Tel * (exp.((E.-μ)./(Constants.kB*Tel)).+1).^2)
end
"""
    dFDdT(Tel::Real, μ::Real, E::Union{Vector{<:Real},Real})
    
    Returns the change in the Fermi distribution with respect to temperature at the given energy value or range.

    # Arguments
    - 'Tel': Thermal electronic temperature
    - 'μ': Chemical Potential
    - 'E': Energy value or range

    # Returns
    - Value or vector of the partial derivative of the Fermi distribution with respect to temperature
"""
@inline function dFDdT(Tel::Real, μ::Real, E::Union{Vector{<:Real},Real})
    return (E.-μ) .* exp.((E.-μ)./(Constants.kB*Tel)) ./ (Constants.kB*Tel^2 * (exp.((E.-μ)./(Constants.kB*Tel)).+1).^2)
end
"""
    dFDdμ(Tel::Real, μ::Real, E::Union{Vector{<:Real},Real})
    
    Returns the change in the Fermi distribution with respect to chemical potential at the given energy value or range.

    # Arguments
    - 'Tel': Thermal electronic temperature
    - 'μ': Chemical Potential
    - 'E': Energy value or range

    # Returns
    - Value or vector of the partial derivative of the Fermi distribution with respect to chemical potential
"""
@inline function dFDdμ(Tel::Real, μ::Real, E::Union{Vector{<:Real},Real})
    return exp.((E.-μ)./(Constants.kB*Tel)) ./ (Constants.kB*Tel * (exp.((E.-μ)./(Constants.kB*Tel)).+1).^2)
end
"""
    build_group_velocity(v_g::Union{Vector{<:Real},Nothing}, FE::Union{Real,Vector{<:Real}}, Conductivity::Bool, conductive_velocity::Symbol, structure::Structure)
    
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
function build_group_velocity(v_g::Union{Vector{<:Real},Nothing}, FE::Union{Real,Vector{<:Real}}, Conductivity::Bool, conductive_velocity::Symbol, structure::Structure)
    if isnothing(v_g)
        if Conductivity == true
            if Conductive_velocity == :fermigas
                if structure.Elemental_System > 1
                    elements = structure.Elemental_System
                    v_g = Vector{Vector{<:Real}}(undef,elements)
                    grids = split_grid(structure.dimension.grid,structure.dimension.InterfaceHeight)
                    for i in 1:Elemental_System
                        length = length(grids[i])
                        v_g[i] = fill(get_fermigas_velocity(Ref(structure.egrid),FE),length)
                    end
                else
                    return get_fermigas_velocity(Ref(structure.egrid),FE)
                end
            elseif Conductive_velocity == :effectiveoneband
                if structure.Elemental_System != 1
                    if structure.Spatial_DOS == true
                        for i in 1:structure.dimension.length
                            v_g[i] = effective_one_band_velocity(structure.DOS[i],egrid,FE[i])
                        end
                    else
                        v_g = fill(zeros(length(structure.egrid)),structure.dimension.length)
                        grids = split_grid(structure.dimension.grid,structure.dimension.InterfaceHeight) 
                        for i in 1:structure.Elemental_System
                            for j in grids[i]
                                v_g[j] = effective_one_band_velocity(structure.DOS[i],egrid,FE[i])
                            end
                        end
                    end
                else
                    if structure.Spatial_DOS == true
                        v_g = effective_one_band_velocity(structure.DOS[end],egrid,FE)
                    else
                        v_g = effective_one_band_velocity(structure.DOS,egrid,FE)
                    end
                end
            end
        else 
            return [0.0]
        end
    else
        if Conductivity == true
            if conductive_velocity == :constant
                if structure.Elemental_System > 1
                    elements = structure.Elemental_System
                    V_g = Vector{Vector{<:Real}}(undef,elements)
                    grids = split_grid(structure.dimension.grid,structure.dimension.InterfaceHeight)
                    for i in 1:Elemental_System
                        length = length(grids[i])
                        V_g[i] = fill(convert_units(v_g),length)
                    end
                else
                    return convert_units(v_g)
                end
            end
        else 
            return [0.0]
        end
    end
end
"""
    get_fermigas_velocity(egrid::Vector{<:Real}, EF::Real)
    
    The analytical free electron gas group velocity, requested by conductive_velocity = :fermigas

    # Arguments
    - 'egrid': Energy grid all distributions are solved on
    - 'FE': The Fermi energy, calculated from get_FermiEnergy

    # Returns
    - The free electron gas group velocity
"""
function get_fermigas_velocity(egrid::Vector{<:Real}, EF::Real)
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
function get_fermigas_dos(egrid::Vector{<:Real}, EF::Real)
    comp1 = 1/(3*π^2)
    comp2 = (2*Constants.me/Constants.ħ^2)^(3/2)
    comp3 = (egrid.+FE).^0.5
    return comp1 .* comp2 .* comp3
end
"""
    effective_one_band_velocity(DOS::spl, egrid::Vector{<:Real}, FE::Real)
    
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
function effective_one_band_velocity(DOS::spl, egrid::Vector{<:Real}, FE::Real)
    k_E = effective_onebandmodel(DOS,egrid,FE)
    v_g = similar(k_E)
    if eltype(v_g) <: AbstractVector
        for j in eachindex(v_g)
            E_k = DataInterpolations.AkimaInterpolation(egrid,k_E,extrapolation = ExtrapolationType.Constant)
            dE_dk = DataInterpolations.derivative.(Ref(E_k),k_E)
            kE_spl = Lightmatter.get_interpolant(egrid,k_E)
            dEdk_spl = Lightmatter.get_interpolant(k_E,dE_dk)
            v_g[j] = dEdk_spl(kE_spl(egrid))
        end
    else
        E_k = DataInterpolations.AkimaInterpolation(egrid,k_E,extrapolation = ExtrapolationType.Constant)
        dE_dk = DataInterpolations.derivative.(Ref(E_k),k_E)
        kE_spl = Lightmatter.get_interpolant(egrid,k_E)
        dEdk_spl = Lightmatter.get_interpolant(k_E,dE_dk)
        v_g = dEdk_spl(kE_spl(egrid))
    end
    return v_g
end
"""
    effective_onebandmodel(DOS::spl, egrid::Vector{<:Real}, FE::Real)
    
    Calculates the dispersion relation within the effective one band model.
    For more details see Mueller & Rethfeld, Phys. Rev. B 87, 035139.

    # Arguments
    - 'DOS': The density-of-states of the system
    - 'egrid': Energy grid all distributions are solved on
    - 'FE': The Fermi energy, calculated from get_FermiEnergy

    # Returns
    - The effective one band model dispersion relation
"""
function effective_onebandmodel(DOS::spl, egrid::Vector{<:Real}, FE::Real)
    if DOS isa Vector{spl}
        k_E=fill(zeros(egrid),length(DOS))
    else
        k_E = zeros(length(egrid))
    end
    
    factor = 3π^2
    int(u,p) = DOS(u)

    if DOS isa Vector{spl}
        for v in eachindex(k_E)
            for (i,E) in enumerate(egrid)
                prob = IntegralProblem(int, -FE, E)
                sol = solve(prob, HCubatureJL(initdiv=100), abstol=1e-8, reltol=1e-8)
                k_E[v][i] = cbrt(factor*sol.u)
            end
        end
    else 
        for (i,E) in enumerate(egrid)
            prob=IntegralProblem(int, -FE, E)
            sol = solve(prob, HCubatureJL(initdiv=100), abstol=1e-8, reltol=1e-8)
            k_E[i] = cbrt(factor*sol.u)
        end
    end
    return k_E
end
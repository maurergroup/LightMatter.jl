"""
    FermiDirac(Tel::Real,μ::Real,E::Union{Vector{<:Real},Real})
    Returns the Fermi distribution at the given temperature and chemical potential for either a single energy
    value or across an energy window - given by whether E is a Real or Vector.
"""
@inline FermiDirac(Tel::Real,μ::Real,E::Union{Vector{<:Real},Real}) = 1 ./(exp.((E.-μ)./(Constants.kB*Tel)).+1)
"""
    dFDdE(Tel::Real,μ::Real,E::Union{Vector{<:Real},Real})
    Returns the derivative of a Fermi distribution with respect to energy at the current energy
    point or across an energy window - given by whether E is a Real or Vector. 
"""
@inline function dFDdE(Tel::Real,μ::Real,E::Union{Vector{<:Real},Real})
    return -exp.((E.-μ)./(Constants.kB*Tel))./(Constants.kB*Tel*(exp.((E.-μ)./(Constants.kB*Tel)).+1).^2)
end
"""
    dFDdT(Tel::Real,μ::Real,E::Union{Vector{<:Real},Real})
    Returns the derivative of a Fermi distribution with respect to temperature at the current energy
    point or across an energy window - given by whether E is a Real or Vector. 
"""
@inline function dFDdT(Tel::Real,μ::Real,E::Union{Vector{<:Real},Real})
    return (E.-μ).*exp.((E.-μ)./(Constants.kB*Tel))./(Constants.kB*Tel^2*(exp.((E.-μ)./(Constants.kB*Tel)).+1).^2)
end
"""
    dFDdμ(Tel::Real,μ::Real,E::Union{Vector{<:Real},Real})
    Returns the derivative of a Fermi distribution with respect to chemical potential at the current energy
    point or across an energy window - given by whether E is a Real or Vector. 
"""
@inline function dFDdμ(Tel::Real,μ::Real,E::Union{Vector{<:Real},Real})
    return exp.((E.-μ)./(Constants.kB*Tel))./(Constants.kB*Tel*(exp.((E.-μ)./(Constants.kB*Tel)).+1).^2)
end

function build_group_velocity(v_g,FE,Conductivity,conductive_velocity,structure)
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
                    v_g = Vector{Vector{<:Real}}(undef,elements)
                    grids = split_grid(structure.dimension.grid,structure.dimension.InterfaceHeight)
                    for i in 1:Elemental_System
                        length = length(grids[i])
                        v_g[i] = fill(convert_units(v_g),length)
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

function get_fermigas_velocity(egrid,EF)
    return sqrt.(2*(egrid.+EF)./Constants.me)
end

function get_fermigas_dos(egrid,FE)
    comp1 = 1/(3*π^2)
    comp2 = (2*Constants.me/Constants.ħ^2)^(3/2)
    comp3 = (egrid.+FE).^0.5
    return comp1.*comp2.*comp3
end

function effective_one_band_velocity(DOS,egrid,FE)
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

function effective_onebandmodel(DOS,egrid,FE)
    if DOS isa Vector{spl}
        k_E=fill(zeros(egrid),length(DOS))
    else
        k_E = zeros(length(egrid))
    end
    
    factor = 3π^2#6*pi #6π^2/σ where σ is a spin factor (2 for electrons)
    int(u,p) = DOS(u)

    if DOS isa Vector{spl}
        for v in eachindex(k_E)
            for (i,E) in enumerate(egrid)
                prob=IntegralProblem(int,-FE,E)
                sol = solve(prob,HCubatureJL(initdiv=100),abstol=1e-8,reltol=1e-8)
                k_E[v][i] = cbrt(factor*sol.u)
            end
        end
    else 
        for (i,E) in enumerate(egrid)
            prob=IntegralProblem(int,-FE,E)
            sol = solve(prob,HCubatureJL(initdiv=100),abstol=1e-8,reltol=1e-8)
            k_E[i] = cbrt(factor*sol.u)
        end
    end
    return k_E
end
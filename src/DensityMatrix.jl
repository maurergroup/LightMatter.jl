#
# This file currently contains all work based on DM Anderson Holstein from NQCDynamics. But want to use some of this as a base
# such as keeping the vonNeumann propagator and making it work with any arbritary Hamilotnian. Porbably only therefore requires
# the final few functions of this
#
function vonNeumann(dρ,ρ,H,t)
    dρ .= 1/(1im * Constants.ħ  )*(H*ρ - ρ*H) #0.6852 is hbar in eVfs
    return nothing
end

function build_dm(N,occ)
    basis = []
    for i in 1:N
        ψ = zeros(N)
        ψ[i]=1
        push!(basis,ψ)
    end

    ρ_states = fill(zeros(N,N),N)
    for i in eachindex(basis)
        ρ_states[i] = occ[i]*basis[i]*basis[i]'
    end
    ρ = sum(ρ_states)
    return Complex.(ρ)
end

function thermal_bath_densitymatrix(H,β,ne)
    thermal_occupation = 1 ./(exp.(diag(H[2:end,2:end])./8.617e-5./β).+1)
    ρ=Matrix(Diagonal(thermal_occupation))
    return Complex.(ρ./tr(ρ).*ne)
end

function build_initaldm(H,β,ne,imp_occ)
    ρ_th = thermal_bath_densitymatrix(H,β,ne)
    ρ = Complex.(zeros(bathstates+1,bathstates+1))
    ρ[1,1] = imp_occ
    ρ[2:end,2:end] .= ρ_th
    return ρ
end

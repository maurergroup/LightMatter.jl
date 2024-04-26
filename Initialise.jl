include("StockEquations.jl")
include("Structs.jl")
using .Equations,.Structs

using DelimitedFiles,Unitful,Integrals,Interpolations,LaTeXStrings
using ForwardDiff,Roots,OrdinaryDiffEq,Plots

function main()
    l,mp,sim,gc,sv=Structs.parameterbuilder()
    Equations.DensityOfStates(mp,sv)
    sv.NumElec=Equations.NumberOfElectrons(0.0,mp.FE,sv,gc)
    TRange=range(0,5000)
    Test=zeros(length(TRange))
    for (i,T) in enumnerate(TRange)
        sv.Tel=T
        sv.μ=ChemicalPotential(sv,gc)
        Test[i]=particleconstant(sv,l)
    end
    display(plot(TRange,Test))
end

main()
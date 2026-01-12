using Pkg
cd(raw"C:\Users\u5522838\.julia\dev\LightMatter")
Pkg.activate(".")
using Unitful, CairoMakie

include("src/Geometry.jl")

println("Testing CoreShell geometry...")
CS = build_CoreShell(25.0, 1.0, 5.0)
println("CoreShell grid size: $(size(CS.grid))")
println("CoreShell mask unique values: $(unique(CS.mask))")
fig_cs = plot_coreshell_circles(CS)
display(fig_cs)

println("\nTesting CoreSatellite geometry...")
CSat = build_CoreSatellite(20.0, 5.0, 1.0)
println("CoreSatellite grid size: $(size(CSat.grid))")
println("CoreSatellite mask unique values: $(unique(CSat.mask))")
fig_csat = plot_coresatellite_structure(CSat)
display(fig_csat)

println("\nTesting AntennaReactor geometry...")
AR = build_AntennaReactor(20.0, 30.0, 1.0)
println("AntennaReactor grid size: $(size(AR.grid))")
println("AntennaReactor mask unique values: $(unique(AR.mask))")
fig_ar = plot_antennareactor_structure(AR)
display(fig_ar)

println("\nAll geometries tested successfully!")

# LightMatter.jl
 
LightMatter.jl is a package dedicated to providing a range of methods for 
simulating the interaction of light with nanoscale materials. The package
provides an ever expanding suite of methods and accomponying approximations
to allow users to build a simulation to their desired specification in terms
of accuracy/parameterization/speed/componets.

The range of methods include,
Classical: Two-Temperature Model
Semi-Classical: E-resolved Boltzmann, Athermal Electron Model
Quantum: Dipole Approximation Single Particle State Hamiltonian propagation

A simple way to get started with LightMatter.jl is to build and simulate a
0D Two-Temperature-Model.
```
using LightMatter, Unitful

“First we define our laser input to the system, this includes laser-material 
parameters such as penetration depth (ϵ), reflectivity (R) etc.”
Las = build_Laser(envelope=:Gaussian, FWHM=50.0, ϕ=50.0, Transport=:optical, ϵ=16.3)

“Next we build our thermal electron and thermal phonon systems with relevant material parameters”
Tel = build_ElectronicTemperature(Enabled=true, Electron_PhononCoupling=true, 
      ElectronicHeatCapacity=:linear, ElectronPhononCouplingValue=:constant, 
      γ = 62.9u"J/m^3/K^2", g = 0.3e17u"J/s/m^3/K")
Tph = build_PhononicTemperature(Enabled=true, Electron_PhononCoupling=true, 
      PhononicHeatCapacity=:constant, Cph = 2.49e6u"J/m^3/K")

“Finally, we combine this into our Simulation”
Sim = build_Simulation(electronictemperature = Tel, phononictemperature = Tph, laser=Las) 
```

The different methods/systems available all have accompanying doc strings on th
arguments/approximations that you can make. 
Coming soon full documentation/tutorials as well as a link to a JuliaCon talk will be
provided to make the package easier to use. 
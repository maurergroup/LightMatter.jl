```@setup logging
@info "Expanding src/getting_started.md..."
start_time = time()
```
# [Getting started](@id getting-started)

To get started we shall assemble a classical two-temperature model(TTM) in 0D for gold. For how 
to extend this to a 1D simulation and include thermal transport see [Two Temperature Model](@ref ttm).
We shall use Unitful to make handling units easier. Unitful quantities can be supplied to any parameter in LightMatter.jl
and the code will convert it to the correct unit system. For more on units see [Unit Management](@ref units)
```@setup gs
using LightMatter, Unitful
```

### Creating our systems

The 3 main systems that comprise a Two-Temperature Model are an electronic temperature,
a phononic temperature and a laser that drives the electronic system. The equation of motion
for the TTM in 0D is given below:
```math
C_\text{el}T_\text{el}(t) = -g(T_\text{el} + T_\text{ph}) + S(t) \\
C_\text{ph}T_\text{ph}(t) = g(T_\text{el} + T_\text{ph})
```
Here, $S(t)$ is the equation for the laser, $C_\text{x}$ is the heatcapacity of the respective 
thermal bath and g is the electron-phonon coupling parameter.

In LightMatter.jl we need to create both thermal baths above as well as define an equation for our
laser. This is all handed by 'build' functions. We shall build the most simple of each of the three
systems in turn. First, the electronic thermal bath.
```@example gs
Tel = build_ElectronicTemperature(Enabled=true, 
      Electron_PhononCoupling=true, ElectronicHeatCapacity=:linear, ElectronPhononCouplingValue=:constant,
      γ = 62.9u"J/m^3/K^2", g = 0.3e17u"J/s/m^3/K")
```
This first system shows us the prototypical setup of building a system in LightMatter.jl. First, we turn the system
on with 'Enabled = true'. Next we decide what approximations we want to use, in this case we are using a constant
electron-phonon coupling value and a linear electronic heat capacity defined as $C_\text{el} = \gamma T_\text{el}$.
Finally we then add the required values as Unitful quantities so we don't have to worry about units. 
For the phonon bath,
```@example gs
Tph = build_PhononicTemperature(Enabled=true, 
      Electron_PhononCoupling=true, PhononicHeatCapacity=:constant, 
      Cph = 2.49e6u"J/m^3/K")
```
The same principle applies here as for the electronic bath, but now we don't need to provide details about the coupling
as the phonon bath assumes that there is an electronic system to couple to and the parameters are located there. We decided
on a constant value of the phonon heat capacity and this is the only value we need to provide.
Last but not least is the laser.
```@example gs
Las = build_Laser(envelope=:Gaussian, FWHM=50.0u"fs", ϕ=10.0u"J/m^2", Transport=:optical, ϵ=16.3e-9u"m")
```
Here we have chosen a gaussian envelope for the laser profile. To see the other options see [Lasers](@ref lasers).
We have provided a full-width half-maximum (FWHM) and a fluence (ϕ). The last two settings describe how the laser
interacts with the material. `:optical` states that the interaction depends on the absorption coefficient (1/ϵ).
ϵ is then the inverse of the absorption coefficient a.k.a the penetration depth.

Once we define all of our systems we combine them into the Simulation.
```@example gs; hide
Sim = build_Simulation(electronictemperature = Tel, phononictemperature = Tph, laser=Las)
```

### Using the Simulation

Now that we have the simulation we can go straight to running it if we desire. But before that we can have a look
at the equations of motion that LightMatter.jl has constructed for our simulation. To do this we run,
```@example gs
eom = LightMatter.function_builder(Sim)
```
This returns us a dictionary for each of the equations of motion. In this case the keys are `"Tel"` and `"Tph"`.
We can have a look at the equation of motion for the electronic thermal bath.
```@example gs
eom["Tel"]
```
This equation looks like what we would expect. We can see the beginning is the equation for a normalised gaussian laser.
This is then followed by the subtraction of the electron-phonon coupling parameter multiplied by the temperature difference
of the thermal baths and this is all divided by the expression for the heat capacity that we expected. 

!!! note

    It's always a good idea when creating a new simulation type to check the different equations of motion and see if they
    are what you expect. Make sure that all of the parameters in the equation have been defined by yourself and saved
    in the correct place in the simulation.

We can now run a simulation by defining a couple initial temperatures and a time span,
```@example gs
initialtemps=Dict("Tel"=>300.0,"Tph"=>300.0)
tspan=(-150.0, 1000.0)
```
Finally, we can now run the simulation. We use the `OrdinaryDiffEq` package for performing the time integration and so
all arguments used within a standard `solve` call can be used here. You can also use the `alg` keyword to define 
any time integration algorithm supported by `OrdinaryDiffEq`. 
```@example gs
sol = run_simulation(Sim, initialtemps, tspan, saveat=1.0, abstol=1e-10, reltol=1e-10)
```
Due to using `RecursiveArrayTools` the solution output can be quite tricky to use. To make this easier for instant
plotting / processing you can call the following to output a dictionary of the solution of each system as well as 
the save times as `Array`'s.
```@example gs
results = LightMatter.seperate_results(sol, initialtemps, Sim)
```
The last thing to complete our first simulation is to save the results to a file. This can be handled by
`post_production`. This will save the results and `Simulation` struct to a HDF5 file as well as can perform
some automated post processing. To see what can be done check out [Outputting](@ref outputs).
```
post_production(sol, "Test.hdf5", initialtemps, [:ThermalFermiDistribution], Sim)
```

### What's next?

Now that we've covered the basics of performing a TTM, we're ready to explore all the capabilities
of LightMatter.jl
All the systems follow the same patterns of enabling, defining approximations and adding parameters
but with a couple of extra notes when you start to perform more complicated simulations. Feel free
to check out the Tutorials or Systems sections to found out more.

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```

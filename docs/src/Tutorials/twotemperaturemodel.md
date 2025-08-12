```@setup logging
@info "Expanding src/Tutorials/twotemperaturemodel.md..."
start_time = time()
```

# [1D Two-Temperature Model Simulation](@id ttm)

In this tutorial we will discuss the theory of the Two-Temperature Model(TTM) as well as how to 
perform a 1D TTM simulation. Due to the computational efficiency of the TTM, it is advised to perform
a 1D simulation where possible as this is the more accurate option. The reason for no 2/3D simulations
in LightMatter.jl is that the laser spot size is often much larger than the nanoscale that this package
is designed to simulate so we can treat all points in this plane as equal. 

## Theory of the TTM

The equations of motion for the 1D TTM are as follows:
```math
\begin{aligned}
C_\text{el}\frac{\partial T_\text{el}(z, t)}{\partial t} &= \nabla(\kappa_\text{el}\nabla T_\text{el}) - g(T_\text{el} + T_\text{ph}) + S(t) \\
C_\text{ph}\frac{\partial T_\text{ph}(z, t)}{\partial t} &= \nabla(\kappa_\text{ph}\nabla T_\text{ph}) + g(T_\text{el} + T_\text{ph})
\end{aligned}
```

The aim of this method is to treat the two systems as thermal baths rather than electronic distributions. This has the advantage
of reducing the computational complexity of propagating a distribution in k or E space into a simple scalar value. However,
this takes the assumption that the systems are always in thermal equilibrium with themselves which isn't the case on the 
short femtosecond time scale after excitation with a laser pulse. To model this non-equilbrium you need to extend your 
level of theory to [AthEM](@ref athem) or the [Boltzmann equation](@ref boltzmann). 

### Thermal Baths

``T_\text{el}`` & ``T_\text{ph}`` are the electronic and phononic thermal baths respectively. 

### Heat Capacities

``C_\text{el}`` & ``C_\text{ph}`` are the heat capacities of the electronic and phononic systems. There is a couple of 
approximations that can be used when modelling both of these in LightMatter.jl. For the electronic heat capacity,
we have a constant (`:constant`), linear (`:linear`) and non-linear (`:nonlinear`) heat capacities, which increase
in accuracy in that order. For a phonon heat capacity we have a constant (`:constant`) and non-linear (`:variable`)
form only. The equations & expressions for all these functions are given below. \

**Electronic Heat Capacities**

| KWARG | Expression/Equation |        
| -----  | ------------------ |
|`:constant` | `:(sim.electronictemperature.γ)` |
|`:linear` | `:(sim.electronictemperature.γ * T_\text{el})` |
|`:nonlinear` | ``\int\frac{\partial f(T\text_{el},μ,ϵ)}{\partial T} \text{DOS}(ϵ)ϵ dϵ`` |
|            | `:(LightMatter.nonlinear_electronheatcapacity(Tel, μ, DOS))` |

Both the constant and linear approximations depend on the γ parameter within the electronic temperature whereas the non-linear
approximation depends on a density-of-states (DOS) and the chemical potential (``μ``). To found out more about DOS in LightMatter.jl
see the [AthEM](@ref athem). μ can be made temperature-dependent by setting `ChemicalPotential` to `true` within `build_Structure`. \

**Phononic Heat Capacities**

| KWARG | Expression/Equation |        
| -----  | ------------------ |
|`:constant` | `:(sim.phononictemperature.Cph)` |
|`:variable` | ``9nk_B(\frac{Tph}{θ})^3 *\int\limits_0^{\theta_D} x^4 \frac{e^x}{(e^x-1)^2} dx`` |
|             | `:(LightMatter.variable_phononheatcapacity(Tph, sim.phononictemperature.n, sim.phononictemperature.θ))` |

The constant approximation depends on the parameter ``Cph`` the user defines when constructing the phononic thermal bath.
In the non-linear case, we use Simpson's rule to determine the heat capacity and this depends on the Debye temperature (``θ``)
and the number of atoms per nm^3 (``n``). The Boltzmann constant is a constant within LightMatter.jl so is handled for you.
To see more check out [Unit Handling](@ref units).

### Electron-Phonon Coupling

``g`` is the electorn-phonon coupling constant in the TTM. This parameter can be either a constant(`:constant`) provided 
during construction or a non-linear parameter (`:variable`) calculated each time-step. \

**Electron-Phonon Coupling Constant**

| KWARG | Expression/Equation |        
| -----  | ------------------ |
|`:constant` | `:(sim.electronictemperature.g)` |
|`:variable` | ``\frac{k_Bλω}{\text{DOS}(μ)ħ}\int \text{DOS}(ϵ)^2 \frac{\partial f(T_\text{el}, μ, ϵ)}{\partial E} dϵ`` |
|            | `:(LightMatter.variable_electronphononcoupling(sim.electronictemperature.λ, sim.electronictemperature.ω, DOS, Tel, μ, Tph, sim.structure.egrid))`|

The new parameters in the `:variable` approximation are the electron-phonon mass enhancement factor (``λ``) and
the second moment of the phonon spectrum (``ω``). To found out more check out [lin2008](@cite).

### Laser

The laser input is defined as S(t) in the equations at the top. It only directly interacts with the electronic system.
To find out about the possible lasers see, [Lasers](@ref lasers).

### Thermal Conductivity

The thermal conductivity in LightMatter.jl is restricted to only the z-axis as currently 2/3D models are not possible.
This means that ``\nabla`` represents ``\frac{\partial}{\partial z}``. The gradients are calculated via finite difference,
primarily central difference except at the edges. The thermal conductivity (``κ``) for electrons and phonons is treated
differently. \
*Electrons* \
``\kappa`` or `sim.electronictemperature.κ` represents the room temperature thermal conductivity of the electrons and 
is scaled temperature dependetly during the simulation via ``κ_\text{final} = \kappa_\text{rt}\frac{T_\text{el}}{T_\text{ph}}``.
The phonon thermal conductivity however is treated as a constant. Currently there are no other options but in the future this
may change. 

---

## Example Simulation

Now that all the theory and approximations of the TTM have been discussed, at least in the context of LightMatter.jl, we
will perform a couple of example simulations to show how to run a 1D simulation, the effect this has on temperatures
compared to a 0D simulation, even at nanofilm thicknesses and also discuss the effect of the different approximations
we can invoke. 

!!! note

    If you want to see a simple 0D Two-Temperature Model (TTM) then please check out [Getting Started](@ref getting-started).


We  are going to use the highest level of approximation for this simulation of a 30nm gold film. Due to the computational 
efficiency of the TTM there is no harm in performing the simulation at this level. In addition, heat capacities and electron
-phonon couplings can be very temperature dependent and unless you know that your material is not in the temperature 
window of your simulation it is best to calculate them explicitly. 

First let's define everything we need outside of the thermal baths.
```
    Dim = build_Dimension(collect(0.0:0.5:30.0))
    Struc = build_Structure(dimension=Dim, chemicalpotential = true, bulk_DOS="DOS/Au_DOS.dat", 
                            bulk_geometry = "DOS/Au_geometry.in")
    Las = build_Laser(envelope=:Gaussian, FWHM=50.0u"fs", ϕ=10.0u"J/m^2", Transport=:optical, ϵ=16.3e-9u"m")
```
Here, we have defined a spatial grid to solve the problem on, this is done by giving a vector to `build_Dimension`.
Next we have provided that struct to the `Structure` struct as well as defined that we want the chemical potential
to update with repsect to the electronic temperature. We have also provided a location to find the DOS and a 
geometry file that was used to calculate the DOS. Finally, we have built a very simple Gaussian laser, the same
as the one used in [Getting Started](@ref getting-started).

!!! note

    All automatic implementations of DOS assembly (interpolation) assume that the DOS and geometry are in
    the FHI-AIMS format. There is an aim to make an automatic conversion in the future but if this is not
    the case then you can supply your own spline fit to the `DOS` KWARG in `build_Structure`. To ensure
    type stability use ``LightMatter.get_interpolant(energy, states)`` to create the spline and make sure
    that the units of your DOS are states/eVnm³.

Next let's build our thermal baths.
```
    Tel = build_ElectronicTemperature(Enabled = true, Electron_PhononCoupling = true, Conductivity = true,
                                      ElectronicHeatCapcity = :nonlinear, ElectronPhononCouplingValue = :variable,
                                      κ = 278u"J/s/m/K", λ = 0.18, ω = 13e-4u"eV^2")
    Tph = build_PhononicTemperature(Enabled = true, Electron_PhononCoupling = true, Conductivity = true,
                                    PhononicHeatCapacity = :variable, κ = 2u"J/s/m/K", θ = 162.3u"K", n = 59.0u"nm^-3")
```
We can now go onto build our `Simulation` by bringing all this together.

```
    Sim = build_Simulation(electronictemperature = Tel, phononictemperature = Tph, laser = Las, structure = Struc)
```

And finally, after defining some initial conditions we can run the simulation
```
    initialtemps=Dict("Tel"=>300.0,"Tph"=>300.0)
    tspan=(-150.0, 1000.0)  
    sol = run_simulation(Sim, initialtemps, tspan, saveat=1.0, abstol=1e-10, reltol=1e-10)
```

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after ``runtime s."
```
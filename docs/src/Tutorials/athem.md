```@setup logging
@info "Expanding src/Tutorials/athem.md..."
start_time = time()
```

# [AthEM Simulation](@id athem)

In this tutorial we will discuss the theory of the Athermal Electron Model (AthEM). To see more 
check out [Uehlein2025](@cite). This model takes the computational efficiency of the TTM but adds a 
athermal electron subsystem to simulate electron-hole pairs. The source of the efficiency compared
to the Boltzmann equation is the use of the relaxation-time approximation without the lack of
energy and particle conservation. This gain in efficiency allows users to do more complex simulations
in terms of dimensionality than the Boltzmann equation would allow, e.g. 1D. The discussion on the
thermal componets of AthEM will be minimal here as it all follows the same as in the [Two-Temperature Model](@ref ttm) 
(TTM) simulation. Here we shall first discuss a 0D AthEM simulation, to see about 1D, check out either
the [Surface DOS](@ref surface-dos) or [Antenna-Reactor](@ref antenna-reactor) tutorials. 

## Theory of AthEM

The athermal electron subsystem of AthEM, is described in a similar way to the Boltzmann equation for
a distribution. In AthEM's case it is written as follows:
```math
    \frac{\partial f^*}{\partial t} = \left.\frac{\partial f^*}{\partial t}\right|_\text{absorb} + \left.\frac{\partial f^*}{\partial t}\right|_\text{e^*e} + \left.\frac{\partial f^*}{\partial t}\right|_\text{e^*ph}
```

Each of the partial differential's on the right-hand side denotes a scattering term with photons, thermal electrons,
and phonons respectively. The athermal electron-athermal electron scattering is set to 0 as in most practical
simulations the number of athermal electrons is so small that they are a minor component. 

The scattering of the athermal electrons with photons is described by the sum of two Fermi's Golden Rule expressions.
One for the hole generation and the other for the electron generation.
```math
    \left.\frac{\partial f^*}{\partial t}\right|_\text{absorb} = \left.\frac{\partial f^*}{\partial t}\right|_\text{h^+} + \left.\frac{\partial f^*}{\partial t}\right|_\text{e^-}  \\
    \left.\frac{\partial f^*}{\partial t}\right|_\text{h^+} = \frac{2\pi V}{\hbar}\text{DOS}(E+hv)|M_{E,E_+}|^2 f(E)[1-f(E+hv)] \\ 
    \left.\frac{\partial f^*}{\partial t}\right|_\text{h^+} = \frac{2\pi V}{\hbar}\text{DOS}(E-hv)|M_{E,E_-}|^2 f(E-hv)[1-f(E)] \\ 
```
Here $V$ is the volume of the cell, everything in LightMatter.jl is per $nm^3$ so this is one in the code's implementation.
$f(E)$ is the current total electronic distribution, $\text{DOS}$ is the electronic density-of-states (DOS), $hv$ is the 
laser photon frequency and |M_{E,E'}|^2 are matrix elements that are set such that the number of electrons and holes
generated are equal as well as that the internal energy of the excitation matches the energy inparted by the laser.
In the future there is plans to add further complexity to the matrix elements. The first thing to now notice, is that
there is an explicit and unescapable dependence on the DOS. This is unlike the TTM which can use approximations
that escape this. Currently there are no keyword arguments that can interact with this component.

```math
\left.\frac{\partial f^*}{\partial t}\right|_\text{e^*e} = - \frac{f^*}{\tau_\text{ee}} + \frac{f^\text{eq} - f^\text{rlx}}{\tau_\text{ee}}
```
Here we have a relaxation time approximation in two parts. Firstly, there is the compoionent that sends the 
athermal system to 0 (the first component) and the second component drives the equilibrium distribution
towards the state with the same internal energy as equilbrium and athermal combined ($f^\text{rlx}$). This
additional effect could be considered as generating the generation of secondary athermal electrons. 

### Thermal Baths

$T_\text{el}$ & $T_\text{ph}$ are the electronic and phononic thermal baths respectively. 

### Heat Capacities

$C_\text{el}$ & $C_\text{ph}$ are the heat capacities of the electronic and phononic systems. There is a couple of 
approximations that can be used when modelling both of these in LightMatter.jl. For the electronic heat capacity,
we have a constant (`:constant`), linear (`:linear`) and non-linear (`:nonlinear`) heat capacities, which increase
in accuracy in that order. For a phonon heat capacity we have a constant (`:constant`) and non-linear (`:variable`)
form only. The equations & expressions for all these functions are given below.
*Electronic Heat Capacities*
`:constant` -> `:(sim.electronictemperature.γ)`
`:linear` -> `:(sim.electronictemperature.γ * T_\text{el})`
`:nonlinear` -> $\int\frac{\partial f(T\text_{el},μ,ϵ)}{\partial T} \text{DOS}(ϵ)ϵ dϵ$ 
            -> `:(LightMatter.nonlinear_electronheatcapacity(Tel, μ, DOS))`

Both the constant and linear approximations depend on the γ parameter within the electronic temperature whereas the non-linear
approximation depends on a density-of-states (DOS) and the chemical potential ($μ$). To found out more about DOS in LightMatter.jl
see the [AthEM](@ref athem). μ can be made temperature-dependent by setting `ChemicalPotential` to `true` within `Structure`.
*Phononic Heat Capacities*
`:constant` -> `:(sim.phononictemperature.Cph)`
`:variable` -> $9nk_B(\frac{Tph}{θ})^3 *\int\limits_0^{\theta_D} x^4 \frac{e^x}{(e^x-1)^2} dx$
             -> `:(LightMatter.variable_phononheatcapacity(Tph, sim.phononictemperature.n, sim.phononictemperature.θ))`
The constant approximation depends on the parameter `Cph` the user defines when constructing the phononic thermal bath.
In the non-linear case, we use Simpson's rule to determine the heat capacity and this depends on the Debye temperature ($θ$)
and the number of atoms per nm^3 ($n$). The Boltzmann constant is a constant within LightMatter.jl so is handled for you.
To see more check out [Unit Handling](@ref units).

### Electron-Phonon Coupling

$g$ is the electorn-phonon coupling constant in the TTM. This parameter can be either a constant(`:constant`) provided 
during construction or a non-linear parameter (`:variable`) calculated each time-step. 
`:constant` -> `:(sim.electronictemperature.g)`
`:variable` -> $\frac{k_Bλω}{\text{DOS}(μ)ħ}\int \text{DOS}(ϵ)^2 \frac{\partial f(T_\text{el}, μ, ϵ)}{\partial E} dϵ
            -> `:(LightMatter.variable_electronphononcoupling(sim.electronictemperature.λ, sim.electronictemperature.ω, DOS, Tel, μ, Tph, sim.structure.egrid))`

The new parameters in the `:variable` approximation are the electron-phonon mass enhancement factor ($λ$) and
the second moment of the phonon spectrum ($ω$). To found out more check out [Lin2008](@cite).

### Laser

The laser input is defined as S(t) in the equations at the top. It only directly interacts with the electronic system.
To find out about the possible lasers see, [Lasers](@ref lasers).

### Thermal Conductivity

The thermal conductivity in LightMatter.jl is restricted to only the z-axis as currently 2/3D models are not possible.
This means that $\nabla$ represents $\frac{\partial}{\partial z}$. The gradients are calculated via finite difference,
primarily central difference except at the edges. The thermal conductivity ($κ$) for electrons and phonons is treated
differently.
*Electrons*
$\kappa$ or `sim.electronictemperature.κ` represents the room temperature thermal conductivity of the electrons and 
is scaled temperature dependetly during the simulation via $κ_\text{final} = \kappa_\text{rt}\frac{T_\text{el}}{T_\text{ph}}$.
The phonon thermal conductivity however is treated as a constant. Currently there are no other options but in the future this
may change. 

##Example Simulation

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
    Struc = build_Structure(dimension=Dim, chemicalpotential = true)
    Las = build_Laser(envelope=:Gaussian, FWHM=50.0u"fs", ϕ=10.0u"J/m^2", Transport=:optical, ϵ=16.3e-9u"m")
```
Here, we have defined a spatial grid to solve the problem on, this is done by giving a vector to `build_Dimension`.
Next we have provided that struct to the `Structure` struct as well as defined that we want the chemical potential
to update with repsect to the electronic temperature. Finally, we have built a very simple Gaussian laser, the same
as the one used in [Getting Started](@ref getting-started).

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
@info "...done after $runtime s."
```
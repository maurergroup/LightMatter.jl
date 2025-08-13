```@setup logging
@info "Expanding src/Tutorials/athem.md..."
start_time = time()
```

# [AthEM Simulation](@id athem)

In this tutorial we will discuss the theory of the Athermal Electron Model (AthEM). To see more 
check out [uehlein2025](@cite). This model takes the computational efficiency of the TTM but adds a 
athermal electron subsystem to simulate electron-hole pairs. The source of the efficiency compared
to the Boltzmann equation is the use of the relaxation-time approximation without the lack of
energy and particle conservation. This gain in efficiency allows users to do more complex simulations
in terms of dimensionality than the Boltzmann equation would allow, e.g. 1D. The discussion on the
thermal componets of AthEM will be minimal here as it all follows the same as in the [Two-Temperature Model](@ref ttm) 
(TTM) simulation. Here we shall first discuss a 0D AthEM simulation, to see about 1D, check out either
the [Surface DOS](@ref surface-dos) or [Antenna-Reactor](@ref antenna-reactor) tutorials. 

## Theory of AthEM

### Athermal Electrons 

The athermal electron subsystem of AthEM, is described in a similar way to the Boltzmann equation for
a distribution. In AthEM's case it is written as follows:
```math
    \frac{\partial f^*(E, t)}{\partial t} = \left.\frac{\partial f^*}{\partial t}\right|_\text{absorb} + \left.\frac{\partial f^*}{\partial t}\right|_\text{el^*el} + \left.\frac{\partial f^*}{\partial t}\right|_\text{el^*ph}
```

Each of the partial differential's on the right-hand side denotes a scattering term with photons, thermal electrons,
and phonons respectively. The athermal electron-athermal electron scattering is set to 0 as in most practical
simulations the number of athermal electrons is so small that they are a minor component. 

The scattering of the athermal electrons with photons is described by the sum of two Fermi's Golden Rule expressions.
One for the hole generation and the other for the electron generation.
```math
\begin{align}
    \left.\frac{\partial f^*}{\partial t}\right|_\text{absorb} &= \left.\frac{\partial f^*}{\partial t}\right|_\text{h^+} + \left.\frac{\partial f^*}{\partial t}\right|_\text{e^-}  \\
    \left.\frac{\partial f^*}{\partial t}\right|_\text{h^+} &= \frac{2\pi V}{\hbar}\text{DOS}(E+hv)|M_{E,E_+}|^2 f(E)[1-f(E+hv)] \\ 
    \left.\frac{\partial f^*}{\partial t}\right|_\text{e^-} &= \frac{2\pi V}{\hbar}\text{DOS}(E-hv)|M_{E,E_-}|^2 f(E-hv)[1-f(E)]
\end{align}
```
Here ``V`` is the volume of the cell, everything in LightMatter.jl is per ``nm^3`` so this is one in the code's implementation.
``f(E)`` is the current total electronic distribution, ``\text{DOS}`` is the electronic density-of-states (DOS), ``hv`` is the 
laser photon frequency and ``|M_{E,E'}|^2`` are matrix elements that are set such that the number of electrons and holes
generated are equal as well as that the internal energy of the excitation matches the energy inparted by the laser.
In the future there is plans to add further complexity to the matrix elements. The first thing to now notice, is that
there is an explicit and unescapable dependence on the DOS. This is unlike the TTM which can use approximations
that escape this. Currently there are no keyword arguments that can interact with this component.

```math
\left.\frac{\partial f^*}{\partial t}\right|_\text{el^*el} = - \frac{f^*}{\tau_\text{ee}} + \frac{f^\text{eq} - f^\text{rlx}}{\tau_\text{ee}}
```
Here we have a relaxation time approximation (RTA) in two parts. Firstly, there is the component that sends the 
athermal system to 0 (the first component) and the second component drives the equilibrium distribution
towards the state with the same internal energy as equilbrium and athermal combined (``f^\text{rlx}``). This
additional effect could be considered as generating the generation of secondary athermal electrons. 

The relaxation time itself has a couple of options for electrons. Either a constant (`:constant`) or Fermi-
Liquid Theory (`:FLT`) relaxation time can be used. The Fermi-Liquid time is considered more acurrate due
to the presence of the energy-dependence on the relaxation time. 
Both are given below \
**Athermal Electron Relaxation Time**

| KWARG | Expression/Equation |        
| -----  | ------------------ |
|`:constant` | `:(sim.athermalelectrons.τ)` |
|`:FLT` | ``\tau\frac{\mu^2}{(E-\mu)^2+(\pi k_B T_\text{el})^2}`` |
|       | `:(sim.athermalelectrons.τ * (μ+sim.athermalelectrons.FE)^2 ./((sim.structure.egrid.-μ).^2 .+ (pi*Constants.kB*Tel)^2))` |

In the FLT there is a new variable, `FE`, which isn't in the equation. This
variable is the Fermi energy and is defined as the difference between the bottom and top of the 
valence band. This can be calculated in LightMatter via `FE = FE_initialization(bulk_DOS)` where 
`bulk_DOS` is the same string you provided to `build_Structure` assuming that the DOS Fermi
energy is at 0.0 eV. This variable corrects for the the fact that in LightMatter.jl we treat the
0 K ``μ`` at 0.0 eV rather than FE. Also, in the case of the FLT, `τ` now is a material-
dependent parameter calculated from, ``τ = 128 / \sqrt{3}\pi\omega_p`` where ``\omega_p`` is the
plasmon-frequency of the material. This you must provide yourselves. 

For the interaction with the phonon system we have,
```math
\left.\frac{\partial f^*}{\partial t}\right|_\text{el^*ph} = - \frac{f^*}{\tau_\text{ep}}
```
This is similar to the electron-electron RTA but the difference is there is no 
driving force on the phonons unlike the thermal electrons. The electron-phonon 
relaxation time can only be treated as (`:constant`). We typically use a relaxation-time
derived from, ``\tau_text{ep} = \tau_\text{fft}hv/k_B \theta_D`` where ``\tau_\text{fft}``
is the free-flight time of electrons and ``\theta_D`` is the Debye temperature.

Now that we have constructed the athermal-electron system let us define how it couples
to both of the thermal baths. 

### Thermal Baths

The equations for how the internal energy of the electronic and phononic thermal systems
are very similar between AthEM and the TTM. In 0D they are 

```math
    \frac{\partial u_\text{el}(t)}{\partial t} = - g(T_\text{el} + T_\text{ph}) + \left.\frac{\partial u_\text{el}}{\partial t}\right|_\text{el-el^*} \\
    \frac{\partial u_\text{ph}(t)}{\partial t} = g(T_\text{el} + T_\text{ph}) + \left.\frac{\partial u_\text{ph}}{\partial t}\right|_\text{ph-el^*} 
```

We have an extra term on the right-hand side which denotes the change in energy from the
relaxation with the athermal electrons. The equation for these extra terms are,
```math
    \left.\frac{\partial u_\text{x}}{\partial t}\right|_\text{x-el^*} = -\int \left.\frac{\partial f^*}{\partial t}\right|_\text{el^*-x} \text{DOS}(E) E dE
```

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after ``runtime s."
```
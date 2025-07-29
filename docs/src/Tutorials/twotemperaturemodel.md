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
C_\text{el}\frac{\partial T_\text{el}(z, t)}{\partial t} = \nabla(\kappa\nabla T_\text{el}) - g(T_\text{el} + T_\text{ph}) + S(t) \\
C_\text{ph}\frac{\partial T_\text{ph}(z, t)}{\partial t} = \nabla(\kappa\nabla T_\text{ph}) + g(T_\text{el} + T_\text{ph})
```

The aim of this method is to treat the two systems as thermal baths rather than electronic distributions. This has the advantage
of reducing the computational complexity of propagating a distribution in k or E space into a simple scalar value. However,
this takes the assumption that the systems are always in thermal equilibrium with themselves which isn't the case on the 
short femtosecond time scale after excitation with a laser pulse. To model this non-equilbrium you need to extend your 
level of theory to [AthEM](@ref athem) or the [Boltzmann equation](@ref boltzmann). 

### Thermal Baths

$T_\text{el}$ & $T_\text{ph}$ are the electronic and phononic thermal baths respectively. 

### Heat Capacities

$C_\text{el}$ & $C_\text{ph}$ are the heat capacities of the electronic and phononic systems. There is a couple of 
approximations that can be used when modelling both of these in LightMatter.jl. For the electronic heat capacity,
we have 


!!! note

    If you want to see a simple Two-Temperature Model (TTM) then please check out [Getting Started](@ref getting-started).

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
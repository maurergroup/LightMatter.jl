```@setup logging
@info "Expanding src/index.md..."
start_time = time()
```
# Introduction

Welcome to the docuemntation for the light-matter interaction package LightMatter.jl.
The documentation will cover all aspectsof using the code including how to setup simulations,
the different systems and approximations that you can choose to implement as well as a tutorial
on how to include your own further approximations into a simulation. We will also try to cover
how to perform all calculations that we have previosuly published and hope that you contribute
to keep the package as useful to the community as possible.

### Objectives of the code

- High performance light-matter simulations at a range of theories
- Provide a toolkit of methods/approximations to build the simulation you desire
- Integrate easily with other great tools found in Julia
- To be an open-source, code-available software for the whole light-matter community.

### Features

Here we provide a list of currently implemented features of the code.
We encourage contributions and extensions of the current methods available

#### Systems 

- [Electronic Thermal Bath](@ref electron-temperature)
- [Phononic Thermal Bath](@ref phonon-temperature)
- [Athermal Electron Subsystem](@ref athermal-electrons)
- [Electronic Distribution](@ref electron-distribution)
- [Phononic Distribution](@ref phonon-distribution)
- [Density Matrix Dynamics](@ref density-matrix)

#### Spatial Features

- 0 or 1D simulations: For a simple example see [Two-Temperature Model](@ref ttm)
- [Spatially resolved density of states](@ref surface-dos)
- [Simulating complex systems](@ref antenna-reactor)

### Dynamics with `DifferentialEquations.jl`

The [`DifferentialEquations`](https://diffeq.sciml.ai/stable/) ecosystem from the
[SciML organisation](https://github.com/SciML/) provides a large library of time integration
algorithms. All dynamics propagation is driven by these algorithms and we provide full support
for all ordinary differential equation methods provided within the ['OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/)
package. 

### Installation

#### 1. Install Julia
Download and install the current stable release from the [Julia website](https://julialang.org/downloads/).
For most platforms `Julia` is provided as a precompiled binary and do not require any installation procedure. However, you need to specify the path to julia or create a symbolic link to the executable that is in your systempath. [Juliaup](https://github.com/JuliaLang/juliaup) is the recommended and easiest
way to manage Julia versions and have a smooth installation experience.

#### 2. Install the package
The package can be installed in the same way as any other registered Julia package:
```julia-repl
pkg> add LightMatter
```

#### 4. Use the package!
```julia-repl
julia> using LightMatter
```
You are now free to proceed to the next section and learn how to use the package.
If you would like you can complete step 5 to double check your installation.

#### 5. Run the tests (optional)

To check the package has been installed correctly and everything is working,
you can execute the tests with:
```julia-repl
pkg> test LightMatter
```

!!! warning

    The tests that ship with this package are currently being developed and
    greatly expanded upon.

### How to use this documentation

The first page to read is the [Getting started](@ref getting-started) section which walks through all the ingredients
needed to perform a classical Two-Tempreature Model simulation. 
After this, the reader is free to explore at their leisure since everything else builds directly
upon sections from the [Getting started](@ref getting-started) page.
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```


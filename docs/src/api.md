```@setup logging
@info "Expanding src/api.md..."
start_time = time()
```

### Simulation Types & Constructors

```@autodocs
Modules = [LightMatter]
Order   = [:type, :function, :constant]
Pages = ["SimulationTypes.jl"]
```

### Athermal Electron Distribution
```@autodocs
Modules = [LightMatter]
Order   = [:type, :function, :constant]
Pages = ["AthermalElectrons.jl"]
```

### Density Matrix

```@autodocs
Modules = [LightMatter]
Order   = [:type, :function, :constant]
Pages = ["DensityMatrix.jl"]
```

### Total Electronic Distribution

```@autodocs
Modules = [LightMatter]
Order   = [:type, :function, :constant]
Pages = ["ElectronicDistribution.jl"]
```

### Thermal Electron Bath

```@autodocs
Modules = [LightMatter]
Order   = [:type, :function, :constant]
Pages = ["ElectronicTemperature.jl"]
```

### Lasers

```@autodocs
Modules = [LightMatter]
Order   = [:type, :function, :constant]
Pages = ["Lasers.jl"]
```

### Phononic Distribution

```@autodocs
Modules = [LightMatter]
Order   = [:type, :function, :constant]
Pages = ["PhononicDistribution.jl"]
```

### Thermal Phonon Bath

```@autodocs
Modules = [LightMatter]
Order   = [:type, :function, :constant]
Pages = ["PhononicTemperature.jl"]
```

### Property Calculation

```@autodocs
Modules = [LightMatter]
Order   = [:type, :function, :constant]
Pages = ["PropertyFunctions.jl"]
```

### Running Dynamics

```@autodocs
Modules = [LightMatter]
Order   = [:type, :function, :constant]
Pages = ["RunDynamics.jl"]
```

### Outputting

```@autodocs
Modules = [LightMatter]
Order   = [:type, :function, :constant]
Pages = ["OutputProcessing.jl"]
```

### DOS / Geometry Functions

```@autodocs
Modules = [LightMatter]
Order   = [:type, :function, :constant]
Pages = ["DOS_Geometry.jl"]
```

#### Unit Handling

```@autodocs
Modules = [LightMatter]
Order   = [:type, :function, :constant]
Pages = ["UnitManagement.jl"]
```

### System Construction

```@autodocs
Modules = [LightMatter]
Order   = [:type, :function, :constant]
Pages = ["SystemConstruction.jl"]
```

### Antenna-Reactor Construction

```@autodocs
Modules = [LightMatter]
Order   = [:type, :function, :constant]
Pages = ["AntennaReactor.jl"]
```

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
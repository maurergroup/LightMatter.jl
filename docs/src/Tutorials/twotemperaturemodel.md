```@setup logging
@info "Expanding src/Tutorials/twotemperaturemodel.md..."
start_time = time()
```

# [1D Two-Temperature Model Simulation](@id ttm)

If you want to see a simple Two-Temperature Model (TTM) then please check out [Getting Started](@ref getting-started).
In this case we shall extend the TTM to be in a single dimension and include thermal transport via a Fick's law
diffusion. 

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
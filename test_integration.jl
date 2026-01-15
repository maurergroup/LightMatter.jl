using LightMatter, SpecialFunctions

# Test 1: x^2
println("Test 1: Integrate x^2 from 0 to 1")
x = collect(0.0:0.001:1.0)
y = x.^2
result = LightMatter.integration_algorithm(y, x)
analytical = 1.0/3.0
error = abs(result - analytical)
println("  Result: $result")
println("  Analytical: $analytical")
println("  Error: $error")
println("  log10(error): $(log10(error))")
println()

# Test 2: sin(x)
println("Test 2: Integrate sin(x) from 0 to π")
x = collect(range(0.0, π, length=513))
y = sin.(x)
result = LightMatter.integration_algorithm(y, x)
analytical = 2.0
error = abs(result - analytical)
println("  Result: $result")
println("  Analytical: $analytical")
println("  Error: $error")
println("  log10(error): $(log10(error))")
println()

# Test 3: exp(-x^2)
println("Test 3: Integrate exp(-x^2) from -3 to 3")
x = collect(range(-3.0, 3.0, length=257))
y = exp.(-x.^2)
result = LightMatter.integration_algorithm(y, x)
analytical = sqrt(π) * erf(3.0)
error = abs(result - analytical)
println("  Result: $result")
println("  Analytical: $analytical")
println("  Error: $error")
println("  log10(error): $(log10(error))")
println()

# Test 4: 1/(1+x^2)
println("Test 4: Integrate 1/(1+x^2) from -1 to 1")
x = collect(-1.0:0.001:1.0)
y = 1.0 ./ (1.0 .+ x.^2)
result = LightMatter.integration_algorithm(y, x)
analytical = π/2
error = abs(result - analytical)
println("  Result: $result")
println("  Analytical: $analytical")
println("  Error: $error")
println("  log10(error): $(log10(error))")

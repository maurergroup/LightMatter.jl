using Documenter
using LightMatter  # ← Replace with your actual package name

makedocs(
    sitename = "LightMatter Documentation",
    format = Documenter.HTML(),
    modules = [LightMatter],
    pages = [
        "Home" => "index.md",
        "Constructing Simulations" => "ConstructingSimulations.md",
        "Running Simulations" => "RunningSimulations.md",
        "Units" => "Units.md",
    ],
    checkdocs=:none
)
using Documenter
using DocumenterCitations
using LightMatter

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"))

@time makedocs(;
    plugins=[bib],
    sitename = "LightMatter.jl",
    modules=[LightMatter, OrdinaryDiffEq, Unitful],
    doctest=false,
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://LightMatter.github.io//stable/",
        #assets=["assets/favicon.ico", "assets/citations.css"],
        size_threshold = 5000*1024,
        example_size_threshold = 8000*1024,
    ),
    authors = "Henry Snowden",
    pages = [
        "Introduction" => "index.md",
        "Getting Started" > "gettingstarted.md",
        "Tutorials" => Any["twotemperaturemodel.md", "athem.md", "boltzmann.md", "densitymatrix.md", "antennareactor.md", "surfaces.md"],
        "Systems" => Any["electronictemperature.md","phononictemperature.md","lasers.md","electronicdistribution.md","phononicdistribution.md",
                         "athermalelectrons.md","densitymatrix.md"],
        "Unit Management" => "Units.md",
        "API" =>"api.md",
        "References" => "references.md"
    ],
    checkdocs=:none
)
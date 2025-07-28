using Documenter
using DocumenterCitations
using LightMatter, OrdinaryDiffEq, Unitful

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"))

function find_all_files(directory)
    files = String[]
    for potential_file in sort(readdir(joinpath(@__DIR__, "src", directory)))
        if potential_file[1] != "." # Ensure we aren't adding any hidden files
            push!(files, joinpath(directory, potential_file))
        end
    end
    return files
end

@time makedocs(;
    plugins=[bib],
    sitename = "LightMatter.jl",
    modules=[LightMatter],
    doctest=false,
    #= format=Documenter.HTML(
        #prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://LightMatter.github.io//stable/",
        #assets=["assets/favicon.ico", "assets/citations.css"],
        #size_threshold = 5000*1024,
        #example_size_threshold = 8000*1024,
    ), =#
    authors = "Henry Snowden",
    pages = [
        "Introduction" => "index.md"
        "Getting Started" => "gettingstarted.md"
        "Tutorials" => find_all_files("Tutorials")
        "Systems" => find_all_files("Systems")
        "Post Processing" => "outputting.md"
        "Unit Management" => "units.md"
        "API" =>"api.md"
        "References" => "references.md"
    ])

if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo="github.com/maurergroup/LightMatter.jl",
        push_preview=true
    )
end
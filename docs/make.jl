push!(LOAD_PATH,"../src/")

using Documenter, Insolation

makedocs(
    sitename="Insolation.jl Documentation",
    doctest = false,
    strict = false,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",),
    modules = [Documenter, Insolation],
    clean = false,
    pages = Any[
        "Home" => "index.md"
    ]
)

deploydocs(
    repo = "github.com/claresinger/Insolation.jl.git",
    target = "build"
)
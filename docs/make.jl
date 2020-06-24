using Documenter, Insolation

makedocs(
    sitename="Insolation.jl",
    doctest = false,
    strict = false,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",),
    modules = [Documenter, Insolation],
    clean = false,
    pages = Any[
        "Home" => "index.md"
        "Zenith Angle Examples" => "ZenithAngleExamples.md"
        "Insolation Examples" => "InsolationExamples.md"
        "Library" => "library.md"
    ]
)

deploydocs(
    repo = "github.com/claresinger/Insolation.jl.git",
    target = "build"
)

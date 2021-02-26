using Documenter, Insolation

# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "100"

makedocs(
    sitename="Insolation.jl",
    doctest = false,
    strict = false,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",),
    modules = [Documenter, Insolation],
    clean = false,
    pages = Any[
        "Home" => "index.md"
        "Zenith Angle Equations" => "ZenithAngleEquations.md"
        "Insolation Examples" => "InsolationExamples.md"
        "APIs" => "library.md"
    ]
)

deploydocs(
    repo = "github.com/claresinger/Insolation.jl.git",
    target = "build"
)

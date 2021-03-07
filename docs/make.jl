using Insolation, Documenter

# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "100"

pages = Any[
        "Home" => "index.md"
        "Zenith Angle Equations" => "ZenithAngleEquations.md"
        "Insolation Examples" => "InsolationExamples.md"
        "APIs" => "library.md"
]

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    collapselevel = 1,
)

makedocs(
    sitename = "Insolation.jl",
    format = format,
    clean = true,
    strict = true,
    modules = [Insolation],
    pages = pages,
)

deploydocs(
    repo = "github.com/claresinger/Insolation.jl.git",
    target = "build",
    push_preview = true,
)
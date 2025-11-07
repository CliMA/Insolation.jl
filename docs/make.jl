using Insolation, Documenter

# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "100"

pages = Any[
    "Home"=>"index.md",
    "Getting Started"=>"GettingStarted.md",
    "Examples"=>"InsolationExamples.md",
    "Milankovitch Cycles"=>"Milankovitch.md",
    "Mathematical Background"=>"SolarGeometry.md",
    "API Reference"=>"library.md",
]

format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true", collapselevel = 1)

makedocs(
    sitename = "Insolation.jl",
    format = format,
    clean = true,
    modules = [Insolation],
    pages = pages,
    checkdocs = :none,
    warnonly = true,  # Don't fail on warnings
)

deploydocs(
    repo = "github.com/CliMA/Insolation.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)

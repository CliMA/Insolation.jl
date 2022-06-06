using Test

push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using Insolation

using Dates
using Statistics
using Roots
using Optim
using DelimitedFiles

using CLIMAParameters
using CLIMAParameters.Planet
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

FT = Float32
@testset "Orbital Params" begin
    include("test_orbit_param.jl")
end
@testset "Types" begin
    include("test_types.jl")
end
@testset "Zenith Angle" begin
    include("test_zenith_angle.jl")
end
@testset "Insolation" begin
    include("test_insolation.jl")
    include("compare_insolation.jl")
end
@testset "Equinox Date" begin
    include("test_perihelion.jl")
    include("test_equinox.jl")
end

FT = Float64
@testset "Orbital Params" begin
    include("test_orbit_param.jl")
end
@testset "Types" begin
    include("test_types.jl")
end
@testset "Zenith Angle" begin
    include("test_zenith_angle.jl")
end
@testset "Insolation" begin
    include("test_insolation.jl")
    include("compare_insolation.jl")
end
@testset "Equinox Date" begin
    include("test_perihelion.jl")
    include("test_equinox.jl")
end

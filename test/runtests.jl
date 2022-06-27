using Test

push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using Insolation

using Dates
using Statistics
using Roots
using Optim
using DelimitedFiles

import CLIMAParameters as CP
import Insolation.Parameters as IP
FT = Float32

include(joinpath(pkgdir(Insolation), "parameters", "create_parameters.jl"))
param_set = create_insolation_parameters(FT)

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
param_set = create_insolation_parameters(FT)

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

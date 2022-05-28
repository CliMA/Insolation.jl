using Test

@testset "Types" begin
    include("test_types.jl")
end

FT = Float64
@testset "Zenith Angle" begin
    include("test_zenith_angle.jl")
end
@testset "Insolation" begin
    include("test_insolation.jl")
end
@testset "Orbits" begin
    include("test_perihelion.jl")
    include("test_equinox.jl")
end

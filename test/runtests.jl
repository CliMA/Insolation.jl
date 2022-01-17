using Test

@testset "Zenith Angle" begin
    include("test_zenith_angle.jl")
end

@testset "Insolation" begin
    include("test_insolation.jl")
end

# @testset "Long-term Variations" begin
#     include("test_equinox.jl")
#     include("test_perihelion.jl")
# end
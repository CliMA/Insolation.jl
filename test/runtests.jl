using Test

@testset "Zenith Angle" begin
    include("test_zenith_angle.jl")
end

@testset "Insolation" begin
    include("test_insolation.jl")
end

# julia --project --color=yes -e 'using Pkg; Pkg.instantiate(); Pkg.build(); Pkg.test();'
using Test

push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Dates
using Statistics
using Roots
using Optim

using Insolation
import Insolation.Parameters as IP
import Insolation.OrbitalDataSplines
import ClimaParams as CP

@testset "Insolation.jl" begin
    for FT_type in (Float32, Float64)
        # Set global variables that test files will access
        global FT = FT_type
        global param_set = IP.InsolationParameters(FT)

        @testset "Type: $FT_type" begin
            @testset "Orbital Params" begin
                include("test_orbit_param.jl")
            end
            @testset "TSIDataSpline" begin
                include("test_tsi_data_spline.jl")
            end
            @testset "Types" begin
                include("test_types.jl")
            end
            @testset "Zenith Angle" begin
                include("test_zenith_angle.jl")
            end
            @testset "Insolation" begin
                include("test_insolation.jl")
            end
            @testset "Daily Insolation" begin
                include("test_daily_insolation.jl")
            end
            @testset "Azimuth Angle" begin
                include("test_azimuth.jl")
            end
            @testset "Physical Consistency" begin
                include("test_physical_consistency.jl")
            end
            @testset "Edge Cases" begin
                include("test_edge_cases.jl")
            end
            @testset "Equinox Date" begin
                include("test_perihelion.jl")
                include("test_equinox.jl")
            end
        end
    end

    # GPU tests (optional - only run if CUDA is available)
    @testset "GPU Compatibility" begin
        include("test_gpu.jl")
    end
end

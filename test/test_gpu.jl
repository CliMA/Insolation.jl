using Test
using Insolation
using Dates
using ClimaParams
using Adapt

# Check if CUDA is available
const CUDA_AVAILABLE = Sys.isunix() && try
    using CUDA
    CUDA.functional()
catch
    false
end

if CUDA_AVAILABLE
    using CUDA

    @testset "GPU Compatibility Tests" begin
        # Test both Float32 and Float64
        for FT in (Float32, Float64)
            @testset "GPU tests with $FT" begin
                # Create parameters on CPU
                params = InsolationParameters(FT)

                # Test data
                date = DateTime(2024, 6, 21, 12, 0, 0)
                lat_cpu = FT(40.0)
                lon_cpu = FT(-105.0)

                # Compute reference on CPU
                F_cpu, S_cpu, μ_cpu, ζ_cpu =
                    insolation(date, lat_cpu, lon_cpu, params)

                # Create orbital data splines
                od_cpu = OrbitalDataSplines()
                od_gpu = adapt(CuArray, od_cpu)

                # Test single value on GPU
                @testset "Single value computation" begin
                    lat_gpu = CuArray([lat_cpu])
                    lon_gpu = CuArray([lon_cpu])

                    # Run on GPU (milankovitch=false, no orbital_data needed)
                    result = insolation.(date, lat_gpu, lon_gpu, params)

                    # Bring back to CPU for comparison
                    F_gpu, S_gpu, μ_gpu, ζ_gpu = Array(result)[1]

                    # Check results match
                    @test F_gpu ≈ F_cpu rtol = 1e-4
                    @test S_gpu ≈ S_cpu rtol = 1e-4
                    @test μ_gpu ≈ μ_cpu rtol = 1e-4
                    @test ζ_gpu ≈ ζ_cpu rtol = 1e-4

                    # Compute with orbital_data and milankovitch = true
                    F_cpu, S_cpu, μ_cpu, ζ_cpu = insolation(
                        date,
                        lat_cpu,
                        lon_cpu,
                        params,
                        od_cpu,
                        milankovitch = true,
                    )
                    result =
                        insolation.(
                            date,
                            lat_gpu,
                            lon_gpu,
                            params,
                            od_gpu,
                            milankovitch = true,
                        )
                    F_gpu, S_gpu, μ_gpu, ζ_gpu = Array(result)[1]

                    @test F_gpu ≈ F_cpu rtol = 1e-4
                    @test S_gpu ≈ S_cpu rtol = 1e-4
                    @test μ_gpu ≈ μ_cpu rtol = 1e-4
                    @test ζ_gpu ≈ ζ_cpu rtol = 1e-4
                end

                # Combinations that we want to try for insolation
                combinations =
                    ((nothing, nothing, false), (od_cpu, od_gpu, true))

                # Test broadcasting over multiple values
                @testset "Broadcasting over multiple locations" begin
                    for (maybe_od_cpu, maybe_od_gpu, milankovitch) in
                        combinations
                        n = 100
                        lats_cpu = FT.(range(-90, 90, length = n))
                        lons_cpu = FT.(range(-180, 180, length = n))

                        # Compute on CPU
                        results_cpu =
                            insolation.(
                                date,
                                lats_cpu,
                                lons_cpu,
                                params,
                                maybe_od_cpu;
                                milankovitch,
                            )

                        # Transfer to GPU
                        lats_gpu = CuArray(lats_cpu)
                        lons_gpu = CuArray(lons_cpu)

                        # Compute on GPU
                        results_gpu =
                            insolation.(
                                date,
                                lats_gpu,
                                lons_gpu,
                                params,
                                maybe_od_gpu;
                                milankovitch,
                            )

                        # Bring back to CPU
                        results_gpu_cpu = Array(results_gpu)

                        # Check all results match
                        for i in 1:n
                            F_cpu, S_cpu, μ_cpu, ζ_cpu = results_cpu[i]
                            F_gpu, S_gpu, μ_gpu, ζ_gpu = results_gpu_cpu[i]

                            @test F_gpu ≈ F_cpu rtol = 1e-4
                            @test S_gpu ≈ S_cpu rtol = 1e-4
                            @test μ_gpu ≈ μ_cpu rtol = 1e-4
                            @test ζ_gpu ≈ ζ_gpu rtol = 1e-4
                        end
                    end
                end

                # Test daily_insolation
                @testset "daily_insolation on GPU" begin
                    for (maybe_od_cpu, maybe_od_gpu, milankovitch) in
                        combinations
                        lat_gpu = CuArray([lat_cpu])

                        # Compute reference on CPU
                        F_daily_cpu, S_daily_cpu, μ_daily_cpu =
                            daily_insolation(
                                date,
                                lat_cpu,
                                params,
                                maybe_od_cpu;
                                milankovitch,
                            )

                        # Run on GPU
                        result =
                            daily_insolation.(
                                date,
                                lat_gpu,
                                params,
                                maybe_od_gpu;
                                milankovitch,
                            )

                        # Bring back to CPU
                        F_daily_gpu, S_daily_gpu, μ_daily_gpu = Array(result)[1]

                        # Check results match
                        @test F_daily_gpu ≈ F_daily_cpu rtol = 1e-5
                        @test S_daily_gpu ≈ S_daily_cpu rtol = 1e-5
                        @test μ_daily_gpu ≈ μ_daily_cpu rtol = 1e-5
                    end
                end

                # Test with different times
                @testset "Multiple dates" begin
                    for (maybe_od_cpu, maybe_od_gpu, milankovitch) in
                        combinations
                        dates = [
                            DateTime(2024, 3, 20, 12, 0, 0),  # Equinox
                            DateTime(2024, 6, 21, 12, 0, 0),  # Summer solstice
                            DateTime(2024, 9, 22, 12, 0, 0),  # Equinox
                            DateTime(2024, 12, 21, 12, 0, 0), # Winter solstice
                        ]
                        lats = FT.([0.0, 45.0, -45.0, 90.0])
                        lons = FT.([0.0, 0.0, 180.0, 0.0])

                        # CPU computation
                        results_cpu =
                            insolation.(
                                dates,
                                lats,
                                lons,
                                params,
                                maybe_od_cpu;
                                milankovitch,
                            )

                        # GPU computation
                        lats_gpu = CuArray(lats)
                        lons_gpu = CuArray(lons)
                        dates_gpu = CuArray(dates)
                        results_gpu =
                            insolation.(
                                dates_gpu,
                                lats_gpu,
                                lons_gpu,
                                params,
                                maybe_od_gpu;
                                milankovitch,
                            )

                        # Compare
                        results_gpu_cpu = Array(results_gpu)
                        for i in 1:length(dates)
                            F_cpu, S_cpu, μ_cpu, ζ_cpu = results_cpu[i]
                            F_gpu, S_gpu, μ_gpu, ζ_gpu = results_gpu_cpu[i]

                            @test F_gpu ≈ F_cpu rtol = 1e-5
                            @test S_gpu ≈ S_cpu rtol = 1e-5
                            @test μ_gpu ≈ μ_cpu rtol = 1e-5
                            @test ζ_gpu ≈ ζ_cpu rtol = 1e-5
                        end
                    end
                end

                # Test edge cases
                @testset "Edge cases on GPU" begin
                    for (maybe_od_cpu, maybe_od_gpu, milankovitch) in
                        combinations
                        # Polar night
                        result_cpu = insolation(
                            DateTime(2024, 12, 21, 12, 0, 0),
                            FT(80.0),
                            FT(0.0),
                            params,
                            maybe_od_cpu;
                            milankovitch,
                        )
                        result_gpu =
                            insolation.(
                                DateTime(2024, 12, 21, 12, 0, 0),
                                CuArray([FT(80.0)]),
                                CuArray([FT(0.0)]),
                                params,
                                maybe_od_gpu;
                                milankovitch,
                            )
                        result_gpu_cpu = Array(result_gpu)[1]

                        @test result_gpu_cpu[1] ≈ result_cpu[1] rtol = 1e-5  # F
                        @test result_gpu_cpu[2] ≈ result_cpu[2] rtol = 1e-5  # S
                        @test result_gpu_cpu[3] ≈ result_cpu[3] rtol = 1e-5  # μ
                    end
                end
            end
        end
    end

    @info "GPU tests completed successfully with CUDA"

else
    @testset "GPU Compatibility Tests" begin
        @test_skip "CUDA not available - skipping GPU tests"
    end
    @info "Skipping GPU tests: CUDA not available or not functional"
end

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
                result_cpu = insolation(date, lat_cpu, lon_cpu, params)
                F_cpu, S_cpu, μ_cpu, ζ_cpu =
                    result_cpu.F, result_cpu.S, result_cpu.μ, result_cpu.ζ

                # Create orbital data splines
                od_cpu = OrbitalDataSplines()
                od_gpu = adapt(CuArray, od_cpu)

                tsi_cpu = TSIDataSpline(FT)
                tsi_gpu = adapt(CuArray, tsi_cpu)

                # Test single value on GPU
                @testset "Single value computation" begin
                    lat_gpu = CuArray([lat_cpu])
                    lon_gpu = CuArray([lon_cpu])

                    # Run on GPU (milankovitch=false, no orbital_data needed)
                    result = insolation.(date, lat_gpu, lon_gpu, params)

                    # Bring back to CPU for comparison
                    result_gpu_arr = Array(result)[1]
                    F_gpu, S_gpu, μ_gpu, ζ_gpu = result_gpu_arr.F,
                    result_gpu_arr.S,
                    result_gpu_arr.μ,
                    result_gpu_arr.ζ

                    # Check results match
                    @test F_gpu ≈ F_cpu rtol = 1e-4
                    @test S_gpu ≈ S_cpu rtol = 1e-4
                    @test μ_gpu ≈ μ_cpu rtol = 1e-4
                    @test ζ_gpu ≈ ζ_cpu rtol = 1e-4

                    # Compute with orbital_data and milankovitch = true
                    milankovitch = true
                    result_cpu2 =
                        insolation(date, lat_cpu, lon_cpu, params, od_cpu, milankovitch)
                    F_cpu, S_cpu, μ_cpu, ζ_cpu =
                        result_cpu2.F, result_cpu2.S, result_cpu2.μ, result_cpu2.ζ
                    result =
                        insolation.(date, lat_gpu, lon_gpu, params, od_gpu, milankovitch)
                    result_gpu_arr = Array(result)[1]
                    F_gpu, S_gpu, μ_gpu, ζ_gpu = result_gpu_arr.F,
                    result_gpu_arr.S,
                    result_gpu_arr.μ,
                    result_gpu_arr.ζ

                    @test F_gpu ≈ F_cpu rtol = 1e-4
                    @test S_gpu ≈ S_cpu rtol = 1e-4
                    @test μ_gpu ≈ μ_cpu rtol = 1e-4
                    @test ζ_gpu ≈ ζ_cpu rtol = 1e-4

                    # Compute with solar variability
                    result_cpu3 = insolation(
                        date,
                        lat_cpu,
                        lon_cpu,
                        params,
                        od_cpu,
                        milankovitch,
                        tsi_cpu,
                    )
                    F_cpu, S_cpu, μ_cpu, ζ_cpu =
                        result_cpu3.F, result_cpu3.S, result_cpu3.μ, result_cpu3.ζ
                    result =
                        insolation.(
                            date,
                            lat_gpu,
                            lon_gpu,
                            params,
                            od_gpu,
                            milankovitch,
                            tsi_gpu,
                        )
                    result_gpu_arr = Array(result)[1]
                    F_gpu, S_gpu, μ_gpu, ζ_gpu = result_gpu_arr.F,
                    result_gpu_arr.S,
                    result_gpu_arr.μ,
                    result_gpu_arr.ζ

                    @test F_gpu ≈ F_cpu rtol = 1e-4
                    @test S_gpu ≈ S_cpu rtol = 1e-4
                    @test μ_gpu ≈ μ_cpu rtol = 1e-4
                    @test ζ_gpu ≈ ζ_cpu rtol = 1e-4
                end

                # Combinations that we want to try for insolation We want to
                # test that the same function on CPU and GPU produces similar
                # results. The function takes the OrbitalDataSplines and
                # TSIDataSpline structs, so we test both functions with the
                # structs on CPU and GPU. Furthermore, there is the milankovitch
                # keyword which we also want to take. The tuple is of the form
                # (od_cpu, od_gpu, milankovitch, tsi_cpu, tsi_gpu) where od and
                # tsi could be nothing.
                combinations = (
                    (nothing, nothing, false, nothing, nothing),
                    (od_cpu, od_gpu, true, nothing, nothing),
                    (od_cpu, od_gpu, true, tsi_cpu, tsi_gpu),
                )

                # Test broadcasting over multiple values
                @testset "Broadcasting over multiple locations" begin
                    for (
                        maybe_od_cpu,
                        maybe_od_gpu,
                        milankovitch,
                        maybe_tsi_cpu,
                        maybe_tsi_gpu,
                    ) in combinations
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
                                maybe_od_cpu,
                                milankovitch,
                                maybe_tsi_cpu,
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
                                maybe_od_gpu,
                                milankovitch,
                                maybe_tsi_gpu,
                            )

                        # Bring back to CPU
                        results_gpu_cpu = Array(results_gpu)

                        # Check all results match
                        for i = 1:n
                            result_cpu_i = results_cpu[i]
                            F_cpu, S_cpu, μ_cpu, ζ_cpu = result_cpu_i.F,
                            result_cpu_i.S,
                            result_cpu_i.μ,
                            result_cpu_i.ζ
                            result_gpu_i = results_gpu_cpu[i]
                            F_gpu, S_gpu, μ_gpu, ζ_gpu = result_gpu_i.F,
                            result_gpu_i.S,
                            result_gpu_i.μ,
                            result_gpu_i.ζ

                            @test F_gpu ≈ F_cpu rtol = 1e-4
                            @test S_gpu ≈ S_cpu rtol = 1e-4
                            @test μ_gpu ≈ μ_cpu rtol = 1e-4
                            @test ζ_gpu ≈ ζ_gpu rtol = 1e-4
                        end
                    end
                end

                # Test daily_insolation
                @testset "daily_insolation on GPU" begin
                    for (
                        maybe_od_cpu,
                        maybe_od_gpu,
                        milankovitch,
                        maybe_tsi_cpu,
                        maybe_tsi_gpu,
                    ) in combinations
                        lat_gpu = CuArray([lat_cpu])

                        # Compute reference on CPU
                        result_daily_cpu = daily_insolation(
                            date,
                            lat_cpu,
                            params,
                            maybe_od_cpu,
                            milankovitch,
                            maybe_tsi_cpu,
                        )
                        F_daily_cpu, S_daily_cpu, μ_daily_cpu =
                            result_daily_cpu.F, result_daily_cpu.S, result_daily_cpu.μ

                        # Run on GPU
                        result =
                            daily_insolation.(
                                date,
                                lat_gpu,
                                params,
                                maybe_od_gpu,
                                milankovitch,
                                maybe_tsi_gpu,
                            )

                        # Bring back to CPU
                        result_daily_gpu = Array(result)[1]
                        F_daily_gpu, S_daily_gpu, μ_daily_gpu =
                            result_daily_gpu.F, result_daily_gpu.S, result_daily_gpu.μ

                        # Check results match
                        @test F_daily_gpu ≈ F_daily_cpu rtol = 1e-5
                        @test S_daily_gpu ≈ S_daily_cpu rtol = 1e-5
                        @test μ_daily_gpu ≈ μ_daily_cpu rtol = 1e-5
                    end
                end

                # Test with different times
                @testset "Multiple dates" begin
                    for (
                        maybe_od_cpu,
                        maybe_od_gpu,
                        milankovitch,
                        maybe_tsi_cpu,
                        maybe_tsi_gpu,
                    ) in combinations
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
                                maybe_od_cpu,
                                milankovitch,
                                maybe_tsi_cpu,
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
                                maybe_od_gpu,
                                milankovitch,
                                maybe_tsi_gpu,
                            )

                        # Compare
                        results_gpu_cpu = Array(results_gpu)
                        for i = 1:length(dates)
                            result_cpu_i = results_cpu[i]
                            F_cpu, S_cpu, μ_cpu, ζ_cpu = result_cpu_i.F,
                            result_cpu_i.S,
                            result_cpu_i.μ,
                            result_cpu_i.ζ
                            result_gpu_i = results_gpu_cpu[i]
                            F_gpu, S_gpu, μ_gpu, ζ_gpu = result_gpu_i.F,
                            result_gpu_i.S,
                            result_gpu_i.μ,
                            result_gpu_i.ζ

                            @test F_gpu ≈ F_cpu rtol = 1e-5
                            @test S_gpu ≈ S_cpu rtol = 1e-5
                            @test μ_gpu ≈ μ_cpu rtol = 1e-5
                            @test ζ_gpu ≈ ζ_cpu rtol = 1e-5
                        end
                    end
                end

                # Test edge cases
                @testset "Edge cases on GPU" begin
                    for (
                        maybe_od_cpu,
                        maybe_od_gpu,
                        milankovitch,
                        maybe_tsi_cpu,
                        maybe_tsi_gpu,
                    ) in combinations
                        # Polar night
                        result_cpu = insolation(
                            DateTime(2024, 12, 21, 12, 0, 0),
                            FT(80.0),
                            FT(0.0),
                            params,
                            maybe_od_cpu,
                            milankovitch,
                            maybe_tsi_cpu,
                        )
                        result_gpu =
                            insolation.(
                                DateTime(2024, 12, 21, 12, 0, 0),
                                CuArray([FT(80.0)]),
                                CuArray([FT(0.0)]),
                                params,
                                maybe_od_gpu,
                                milankovitch,
                                maybe_tsi_gpu,
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

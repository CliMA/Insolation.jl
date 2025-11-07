date = Dates.DateTime(2020, 2, 20, 11, 11, 0)
lon, lat = [FT(80.0), FT(20.0)]

od = Insolation.OrbitalDataSplines()

# Test with fixed epoch parameters (no eot correction)
F, S, μ, ζ = insolation(date, lat, lon, param_set; eot_correction = false)
@test typeof(F) == FT
@test typeof(S) == FT
@test typeof(μ) == FT
@test typeof(ζ) == FT

# Test with Milankovitch cycles (no eot correction)
F, S, μ, ζ =
    insolation(date, lat, lon, param_set, od; milankovitch = true, eot_correction = false)
@test typeof(F) == FT
@test typeof(S) == FT
@test typeof(μ) == FT
@test typeof(ζ) == FT

# Test with fixed epoch parameters (with eot correction)
F, S, μ, ζ = insolation(date, lat, lon, param_set; eot_correction = true)
@test typeof(F) == FT
@test typeof(S) == FT
@test typeof(μ) == FT
@test typeof(ζ) == FT

# Test with Milankovitch cycles (with eot correction)
F, S, μ, ζ =
    insolation(date, lat, lon, param_set, od; milankovitch = true, eot_correction = true)
@test typeof(F) == FT
@test typeof(S) == FT
@test typeof(μ) == FT
@test typeof(ζ) == FT

# Test Base.broadcastable for InsolationParameters
# This ensures InsolationParameters behaves as a scalar in broadcasting
@testset "InsolationParameters Broadcasting" begin
    # Test that broadcastable returns a tuple (makes it a scalar)
    @test Base.broadcastable(param_set) isa Tuple
    @test length(Base.broadcastable(param_set)) == 1
    @test Base.broadcastable(param_set)[1] === param_set

    # Test broadcasting without Ref() - params should act as scalar
    dates = [Dates.DateTime(2024, 6, 21), Dates.DateTime(2024, 12, 21)]
    lats = FT.([40.0, 60.0])
    lons = FT.([-105.0, 25.0])

    # This should work without Ref() because of Base.broadcastable
    results = insolation.(dates, lats, lons, param_set)

    # Verify we got 2 results (one per date/lat/lon combo)
    @test length(results) == 2

    # Verify each result is a tuple of (F, S, μ, ζ)
    for result in results
        @test result isa Tuple
        @test length(result) == 4
        F_test, S_test, μ_test, ζ_test = result
        @test F_test isa FT
        @test S_test isa FT
        @test μ_test isa FT
        @test ζ_test isa FT
    end

    # Test that results match single calls
    F1, S1, μ1, ζ1 = insolation(dates[1], lats[1], lons[1], param_set)
    @test results[1][1] == F1
    @test results[1][2] == S1
    @test results[1][3] == μ1
    @test results[1][4] == ζ1
end

## Test edge cases and boundary conditions. 
## Tests for Earth's current orbital parameters
rtol = 1e-2
atol = 1e-4

od = Insolation.OrbitalDataSplines()

@testset "Edge Cases - Exact Poles" begin
    # Test at exact north and south poles
    date_summer = Dates.DateTime(2000, 6, 21, 12, 0, 0)
    date_winter = Dates.DateTime(2000, 12, 21, 12, 0, 0)
    lon = FT(0.0)

    # North Pole
    lat_np = FT(90.0)

    # Should not error
    (; F, S, μ, ζ) = insolation(date_summer, lat_np, lon, param_set)
    F_np_s, S_np_s, μ_np_s, ζ_np_s = F, S, μ, ζ
    (; F, S, μ, ζ) = insolation(date_winter, lat_np, lon, param_set)
    F_np_w, S_np_w, μ_np_w, ζ_np_w = F, S, μ, ζ

    @test all(isfinite.([F_np_s, S_np_s, μ_np_s, ζ_np_s]))
    @test all(isfinite.([F_np_w, S_np_w, μ_np_w, ζ_np_w]))

    # North pole should have more insolation in NH summer
    @test F_np_s > F_np_w
    @test F_np_w ≈ 0 atol = atol

    # South Pole
    lat_sp = FT(-90.0)

    (; F, S, μ, ζ) = insolation(date_summer, lat_sp, lon, param_set)
    F_sp_s, S_sp_s, μ_sp_s, ζ_sp_s = F, S, μ, ζ
    (; F, S, μ, ζ) = insolation(date_winter, lat_sp, lon, param_set)
    F_sp_w, S_sp_w, μ_sp_w, ζ_sp_w = F, S, μ, ζ

    @test all(isfinite.([F_sp_s, S_sp_s, μ_sp_s, ζ_sp_s]))
    @test all(isfinite.([F_sp_w, S_sp_w, μ_sp_w, ζ_sp_w]))

    # South pole should have more insolation in NH winter (SH summer)
    @test F_sp_w > F_sp_s
    @test F_sp_s ≈ 0 atol = atol
end

@testset "Edge Cases - Date Line" begin
    # Test at international date line (180° longitude)
    date = Dates.DateTime(2000, 6, 21, 12, 0, 0)
    lat = FT(0.0)

    lon_180 = FT(180.0)
    lon_neg180 = FT(-180.0)

    # Both should give same results (same location)
    (; F, S, μ, ζ) = insolation(date, lat, lon_180, param_set)
    F1, S1, μ1, ζ1 = F, S, μ, ζ
    (; F, S, μ, ζ) = insolation(date, lat, lon_neg180, param_set)
    F2, S2, μ2, ζ2 = F, S, μ, ζ

    @test F1 ≈ F2 rtol = rtol
    @test S1 ≈ S2 rtol = rtol
    @test μ1 ≈ μ2 rtol = rtol
    # Azimuth might differ by 2π due to wrapping, check modulo
    @test mod(ζ1, 2π) ≈ mod(ζ2, 2π) rtol = rtol
end

@testset "Edge Cases - Extreme Dates" begin
    # Test with dates far from epoch (but within spline range)
    lat = FT(45.0)
    lon = FT(0.0)

    # Far past
    date_past = Dates.DateTime(1900, 6, 21, 12, 0, 0)
    milankovitch = true
    (; F, S, μ, ζ) = insolation(date_past, lat, lon, param_set, od, milankovitch)
    F_past, S_past, μ_past, ζ_past = F, S, μ, ζ

    @test all(isfinite.([F_past, S_past, μ_past, ζ_past]))
    @test F_past > 0
    @test S_past > 0
    @test 0 <= μ_past <= 1

    # Far future
    date_future = Dates.DateTime(2100, 6, 21, 12, 0, 0)
    milankovitch = true
    (; F, S, μ, ζ) = insolation(date_future, lat, lon, param_set, od, milankovitch)
    F_future, S_future, μ_future, ζ_future = F, S, μ, ζ

    @test all(isfinite.([F_future, S_future, μ_future, ζ_future]))
    @test F_future > 0
    @test S_future > 0
    @test 0 <= μ_future <= 1
end

@testset "Edge Cases - Midnight" begin
    # Test at exactly midnight (potential edge case for hour angle)
    date = Dates.DateTime(2000, 6, 21, 0, 0, 0)
    lat = FT(45.0)
    lon = FT(0.0)

    (; F, S, μ, ζ) = insolation(date, lat, lon, param_set)

    # Should work without errors
    @test all(isfinite.([F, S, μ, ζ]))

    # At midnight, insolation should be zero
    @test F ≈ 0 atol = atol
end

@testset "Edge Cases - Equinoxes" begin
    # Test at exact equinox dates
    equinox_spring = Dates.DateTime(2000, 3, 20, 12, 0, 0)
    equinox_fall = Dates.DateTime(2000, 9, 22, 12, 0, 0)

    lat = FT(0.0)  # Equator
    lon = FT(0.0)

    (; F, μ) = insolation(equinox_spring, lat, lon, param_set)
    F_spring, μ_spring = F, μ
    (; F, μ) = insolation(equinox_fall, lat, lon, param_set)
    F_fall, μ_fall = F, μ

    # At equinoxes, equator should have sun nearly overhead at local noon 
    # (exactly overhead at solar noon, but this is not exactly what we are testing here)
    @test μ_spring > 0.95
    @test μ_fall > 0.95

    # Insolations should be similar (within 5% due to Earth-Sun distance)
    @test F_spring ≈ F_fall rtol = 0.05
end

@testset "Edge Cases - Solstices" begin
    # Test at exact solstice dates
    solstice_summer = Dates.DateTime(2000, 6, 21, 12, 0, 0)
    solstice_winter = Dates.DateTime(2000, 12, 21, 12, 0, 0)

    # At Tropic of Cancer
    lat_toc = FT(rad2deg(IP.obliq_epoch(param_set)))
    lon = FT(0.0)

    (; F, μ) = insolation(solstice_summer, lat_toc, lon, param_set)
    F_summer, μ_summer = F, μ

    # Sun should be nearly overhead at summer solstice on Tropic of Cancer
    @test μ_summer > 0.98

    # At Tropic of Capricorn
    lat_toc_s = FT(-rad2deg(IP.obliq_epoch(param_set)))
    (; F, μ) = insolation(solstice_winter, lat_toc_s, lon, param_set)
    F_winter, μ_winter = F, μ

    # Sun should be nearly overhead at winter solstice on Tropic of Capricorn
    @test μ_winter > 0.98
end

@testset "Edge Cases - Prime Meridian" begin
    # Test at exact 0° longitude
    date = Dates.DateTime(2000, 6, 21, 12, 0, 0)
    lat = FT(45.0)
    lon = FT(0.0)

    (; F, S, μ, ζ) = insolation(date, lat, lon, param_set)

    @test all(isfinite.([F, S, μ, ζ]))
    @test F > 0

    # At solar noon on prime meridian, azimuth should be southerly
    @test abs(ζ - 3π / 2) < 0.2  # Within ~11 degrees of south
end

@testset "Edge Cases - Leap Day" begin
    # Test on February 29 (leap day)
    date_leap = Dates.DateTime(2000, 2, 29, 12, 0, 0)
    lat = FT(45.0)
    lon = FT(0.0)

    (; F, S, μ, ζ) = insolation(date_leap, lat, lon, param_set)

    # Should handle leap day without issues
    @test all(isfinite.([F, S, μ, ζ]))
    @test F > 0

    # Should be similar to day before and after
    date_before = Dates.DateTime(2000, 2, 28, 12, 0, 0)
    date_after = Dates.DateTime(2000, 3, 1, 12, 0, 0)

    (; F) = insolation(date_before, lat, lon, param_set)
    F_before = F
    (; F) = insolation(date_after, lat, lon, param_set)
    F_after = F

    # Should be reasonably close (within 10% due to changing solar position)
    @test 0.9 * F_before < F < 1.1 * F_after
end

@testset "Edge Cases - Daily vs Instantaneous" begin
    # For polar night/day, compare daily and instantaneous results

    # North pole polar night
    date = Dates.DateTime(2000, 12, 21)
    lat = FT(90.0)
    lon = FT(0.0)

    # Daily should be near zero
    (; F, μ) = daily_insolation(date, lat, param_set)
    F_daily, μ_daily = F, μ
    @test F_daily ≈ 0 atol = atol
    @test μ_daily ≈ 0 atol = atol

    # Instantaneous at noon should also be near zero
    date_noon = date + Dates.Hour(12)
    (; F, μ) = insolation(date_noon, lat, lon, param_set)
    F_inst, μ_inst = F, μ
    @test F_inst ≈ 0 atol = atol
    @test μ_inst ≈ 0 atol = atol
end

@testset "Edge Cases - All Longitudes at Same Time" begin
    # Test that all longitudes give valid results at same UTC time
    date = Dates.DateTime(2000, 6, 21, 12, 0, 0)
    lat = FT(45.0)

    longitudes = FT.(range(-180, stop = 180, length = 37))

    for lon in longitudes
        (; F, S, μ, ζ) = insolation(date, lat, lon, param_set)

        # All should give valid finite results
        @test all(isfinite.([F, S, μ, ζ]))
        @test F >= 0
        @test S > 0
        @test 0 <= μ <= 1
        @test 0 <= ζ <= 2π
    end
end

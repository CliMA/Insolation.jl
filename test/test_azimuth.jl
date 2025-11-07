## Test azimuth angle calculations
rtol = 1e-2
atol = 1e-3

od = Insolation.OrbitalDataSplines()

@testset "Azimuth - Solar Noon" begin
    # At solar noon, sun should be due south in Northern Hemisphere
    # Azimuth convention: ζ = 0 is due East, increasing counter-clockwise
    # So due South is ζ = 3π/2

    date = Dates.DateTime(2000, 6, 21, 12, 0, 0)
    lon = FT(0.0)
    lat = FT(45.0)  # Northern hemisphere

    F, S, μ, ζ = insolation(date, lat, lon, param_set)

    # At solar noon, azimuth should be approximately 3π/2 (due south)
    @test ζ ≈ 3π / 2 rtol = rtol

    # Test in Southern Hemisphere
    lat_s = FT(-45.0)
    F_s, S_s, μ_s, ζ_s = insolation(date, lat_s, lon, param_set)

    # In SH, sun at noon should be due north (π/2)
    @test ζ_s ≈ π / 2 rtol = rtol
end

@testset "Azimuth - Sunrise and Sunset" begin
    # Test azimuth at sunrise and sunset
    date = Dates.DateTime(2000, 3, 20)  # Near equinox
    lon = FT(0.0)
    lat = FT(0.0)  # Equator

    # Sunrise (approximately 6 AM)
    sunrise = Dates.DateTime(2000, 3, 20, 6, 0, 0)
    F_sr, S_sr, μ_sr, ζ_sr = insolation(sunrise, lat, lon, param_set)

    # At equinox on equator, sunrise should be due east (ζ ≈ 0 or 2π)
    ζ_sr_normalized = mod(ζ_sr, 2π)
    # Check if close to 0 or close to 2π
    east_test = (abs(ζ_sr_normalized) < 1e-2) || (abs(ζ_sr_normalized - 2π) < 1e-2)
    @test east_test

    # Sunset (approximately 6 PM)
    sunset = Dates.DateTime(2000, 3, 20, 18, 0, 0)
    F_ss, S_ss, μ_ss, ζ_ss = insolation(sunset, lat, lon, param_set)

    # At equinox on equator, sunset should be due west (ζ ≈ π)
    @test ζ_ss ≈ π atol = 1e-2
end

@testset "Azimuth - Bounds Check" begin
    # Azimuth should always be between 0 and 2π
    dates = [Dates.DateTime(2000, 1, 1, h, 0, 0) for h = 0:23]

    lat = FT(45.0)
    lon = FT(0.0)

    for date in dates
        _, _, _, ζ = insolation(date, lat, lon, param_set)
        @test 0 <= ζ <= 2π
    end
end

@testset "Azimuth - Morning vs Afternoon" begin
    # Test azimuth progression through the day
    # Azimuth convention: 0 = East, increasing counter-clockwise
    # So: East (0) → North (π/2) → West (π) → South (3π/2) → East (2π)
    # But sun moves: East (morning) → South (noon) → West (afternoon)
    # This means azimuth: ~0 or 2π (morning) → 3π/2 (noon) → π (afternoon)
    # So azimuth DECREASES through the day in NH!

    date_base = Dates.DateTime(2000, 6, 21)
    lat = FT(45.0)
    lon = FT(0.0)

    # Morning
    morning = date_base + Dates.Hour(8)
    _, _, _, ζ_morning = insolation(morning, lat, lon, param_set)

    # Noon
    noon = date_base + Dates.Hour(12)
    _, _, _, ζ_noon = insolation(noon, lat, lon, param_set)

    # Afternoon
    afternoon = date_base + Dates.Hour(16)
    _, _, _, ζ_afternoon = insolation(afternoon, lat, lon, param_set)

    # Morning should be in eastern quadrant (near 0 or 2π)
    @test abs(ζ_morning - 2π) < 1e-1 || abs(ζ_morning) < 1e-1

    # Noon should be southerly (near 3π/2)
    @test abs(ζ_noon - 3π / 2) < 1e-1

    # Afternoon should be westerly (near π)
    @test abs(ζ_afternoon - π) < 1e-1
end

@testset "Azimuth - Different Longitudes" begin
    date = Dates.DateTime(2000, 6, 21, 12, 0, 0)  # Noon at prime meridian
    lat = FT(45.0)

    # Prime meridian
    lon1 = FT(0.0)
    _, _, _, ζ1 = insolation(date, lat, lon1, param_set)

    # 30° East - solar noon was earlier (should be past noon)
    lon2 = FT(30.0)
    _, _, _, ζ2 = insolation(date, lat, lon2, param_set)

    # 30° West - solar noon is later (should be before noon)
    lon3 = FT(-30.0)
    _, _, _, ζ3 = insolation(date, lat, lon3, param_set)

    # Should be close to south at prime meridian
    @test ζ1 ≈ 3π / 2 rtol = 0.05

    # All should be in reasonable range
    @test 0 <= ζ1 <= 2π
    @test 0 <= ζ2 <= 2π
    @test 0 <= ζ3 <= 2π
end

@testset "Azimuth - High Latitudes" begin
    # Test azimuth behavior at high latitudes
    date = Dates.DateTime(2000, 6, 21, 12, 0, 0)
    lon = FT(0.0)

    # Arctic
    lat_arctic = FT(70.0)
    _, _, _, ζ_arctic = insolation(date, lat_arctic, lon, param_set)

    @test 0 <= ζ_arctic <= 2π

    # Antarctic  
    lat_antarctic = FT(-70.0)
    _, _, _, ζ_antarctic = insolation(date, lat_antarctic, lon, param_set)

    @test 0 <= ζ_antarctic <= 2π
end

@testset "Azimuth - Type Stability" begin
    # Test that azimuth has correct type for both Float32 and Float64
    date = Dates.DateTime(2000, 6, 21, 12, 0, 0)
    lon = FT(0.0)
    lat = FT(45.0)

    _, _, _, ζ = insolation(date, lat, lon, param_set)

    @test typeof(ζ) == FT
end

@testset "Azimuth and Zenith Angles - NOAA Reference Values" begin
    # Test against NOAA Solar Calculator reference values
    # Reference: https://gml.noaa.gov/grad/solcalc/azel.html
    # Tolerances: 1 degree on azimuth, 5e-3 on cos(SZA)

    # Test case 1: Northern hemisphere, autumn
    date1 = Dates.DateTime(2025, 10, 22, 13, 30, 0)
    lat1 = FT(40.0)
    lon1 = FT(-15.0)
    F1, S1, μ1, ζ1 = insolation(date1, lat1, lon1, param_set)

    # Convert azimuth from radians to degrees for comparison
    ζ1_deg = rad2deg(ζ1)
    @test abs(ζ1_deg - 255.83) < 1.0  # 1 degree tolerance
    @test abs(μ1 - 0.6112) < 5e-3     # 5e-3 tolerance on cos(SZA)

    # Test case 2: Southern hemisphere, spring
    date2 = Dates.DateTime(2025, 10, 22, 13, 30, 0)
    lat2 = FT(-50.0)
    lon2 = FT(-15.0)
    F2, S2, μ2, ζ2 = insolation(date2, lat2, lon2, param_set)

    ζ2_deg = rad2deg(ζ2)
    @test abs(ζ2_deg - 107.61) < 1.0
    @test abs(μ2 - 0.7677) < 5e-3

    # Test case 3: Southern hemisphere, autumn
    date3 = Dates.DateTime(2025, 5, 22, 10, 30, 0)
    lat3 = FT(-50.0)
    lon3 = FT(-15.0)
    F3, S3, μ3, ζ3 = insolation(date3, lat3, lon3, param_set)

    ζ3_deg = rad2deg(ζ3)
    @test abs(ζ3_deg - 55.05) < 1.0
    @test abs(μ3 - 0.2161) < 5e-3

    # Test case 4: Northern hemisphere, spring
    date4 = Dates.DateTime(2025, 5, 22, 23, 30, 0)
    lat4 = FT(43.0)
    lon4 = FT(105.0)
    F4, S4, μ4, ζ4 = insolation(date4, lat4, lon4, param_set)

    ζ4_deg = rad2deg(ζ4)
    @test abs(ζ4_deg - 10.99) < 1.0
    @test abs(μ4 - 0.3394) < 5e-3
end

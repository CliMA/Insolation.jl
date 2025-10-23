## Test daily_insolation API
rtol = 1e-2
atol = 1e-4

od = Insolation.OrbitalDataSplines()

@testset "Daily Insolation - Basic Functionality" begin
    # Test with Milankovitch cycles
    date = Dates.DateTime(2000, 6, 21, 12, 0, 0)  # Summer solstice
    lat = FT(45.0)

    F, S, μ = daily_insolation(date, lat, param_set, od; milankovitch = true)

    # Check return types
    @test typeof(F) == FT
    @test typeof(S) == FT
    @test typeof(μ) == FT

    # Check physical bounds
    @test F >= 0
    @test S >= 0
    @test 0 <= μ <= 1

    # Test without Milankovitch cycles (epoch parameters)
    F_epoch, S_epoch, μ_epoch =
        daily_insolation(date, lat, param_set; milankovitch = false)

    @test typeof(F_epoch) == FT
    @test typeof(S_epoch) == FT
    @test typeof(μ_epoch) == FT

    # For year 2000 (near epoch), both methods should give similar results
    @test F ≈ F_epoch rtol = 0.02
    @test S ≈ S_epoch rtol = 0.02
    @test μ ≈ μ_epoch rtol = 0.02
end

@testset "Daily Insolation - Seasonal Variation" begin
    lat = FT(45.0)  # Northern hemisphere

    # Summer solstice - higher insolation in NH
    summer = Dates.DateTime(2000, 6, 21)
    F_summer, _, _ = daily_insolation(summer, lat, param_set)

    # Winter solstice - lower insolation in NH
    winter = Dates.DateTime(2000, 12, 21)
    F_winter, _, _ = daily_insolation(winter, lat, param_set)

    # Test equinoxes - should be between summer and winter
    spring = Dates.DateTime(2000, 3, 20)
    F_spring, _, _ = daily_insolation(spring, lat, param_set)

    # Summer should have more insolation than winter in NH
    @test F_winter < F_spring < F_summer
end

@testset "Daily Insolation - Hemispheric Symmetry" begin
    date = Dates.DateTime(2000, 6, 21)  # NH summer solstice
    lat_north = FT(45.0)
    lat_south = FT(-45.0)

    F_north, _, _ = daily_insolation(date, lat_north, param_set)
    F_south, _, _ = daily_insolation(date, lat_south, param_set)

    # NH should have more insolation than SH at NH summer solstice
    @test F_north > F_south

    # Test at equinox - should be symmetric
    equinox = Dates.DateTime(2000, 3, 20)
    F_north_eq, _, _ = daily_insolation(equinox, lat_north, param_set)
    F_south_eq, _, _ = daily_insolation(equinox, lat_south, param_set)

    @test F_north_eq ≈ F_south_eq rtol = rtol
end

@testset "Daily Insolation - Polar Regions" begin
    # North pole during polar day (summer) - use latitude below pole to avoid singularity
    date_summer = Dates.DateTime(2000, 6, 21)
    lat_np = FT(85.0)  # Near pole but not exactly at it
    F_np_summer, _, μ_np_summer =
        daily_insolation(date_summer, lat_np, param_set)

    # Should have positive insolation during polar day
    @test F_np_summer > 0
    @test μ_np_summer > 0

    # North pole during polar night (winter)
    date_winter = Dates.DateTime(2000, 12, 21)
    F_np_winter, _, μ_np_winter =
        daily_insolation(date_winter, lat_np, param_set)

    # Should have zero or near-zero insolation during polar night
    @test F_np_winter ≈ 0 atol = 1.0  # Small numerical tolerance
    @test μ_np_winter ≈ 0 atol = 0.01

    # South pole - opposite pattern
    lat_sp = FT(-85.0)  # Near pole but not exactly at it
    F_sp_summer, _, μ_sp_summer =
        daily_insolation(date_summer, lat_sp, param_set)
    F_sp_winter, _, μ_sp_winter =
        daily_insolation(date_winter, lat_sp, param_set)

    # SH winter when NH summer
    @test F_sp_summer ≈ 0 atol = 1.0
    @test μ_sp_summer ≈ 0 atol = 0.01
    # SH summer when NH winter  
    @test F_sp_winter > 0
    @test μ_sp_winter > 0
end

@testset "Daily Insolation - Equator" begin
    lat_eq = FT(0.0)

    # Test throughout the year
    dates = [
        Dates.DateTime(2000, 1, 1),
        Dates.DateTime(2000, 4, 1),
        Dates.DateTime(2000, 7, 1),
        Dates.DateTime(2000, 10, 1),
    ]

    insolations = [daily_insolation(d, lat_eq, param_set)[1] for d in dates]

    # All should be positive
    @test all(insolations .> 0)

    # Equatorial insolation should be relatively constant
    # but with some variation due to Earth-Sun distance
    @test maximum(insolations) / minimum(insolations) < 1.15
end

@testset "Daily Insolation - Annual Mean" begin
    # Test that annual mean insolation ≈ TSI/4
    nlats = 37
    ndays = 73

    lats = FT.(collect(range(-90, stop = 90, length = nlats)))
    days = Array{Int}(round.(collect(range(0, stop = 365, length = ndays))))

    F_arr = zeros(FT, ndays, nlats)

    for (i, day) in enumerate(days)
        for (j, lat) in enumerate(lats)
            date = Dates.DateTime(2000, 1, 1) + Dates.Day(day)
            F_arr[i, j], _, _ = daily_insolation(date, lat, param_set)
        end
    end

    # Compute area-weighted global mean
    temporal_mean = mean(F_arr, dims = 1)
    area_weights = abs.(cosd.(lats))
    global_mean = sum(temporal_mean .* area_weights') / sum(area_weights)

    # Should equal TSI/4
    @test global_mean ≈ IP.tot_solar_irrad(param_set) / 4 rtol = rtol
end

@testset "Daily Insolation - Consistency with Base Function" begin
    # Test that high-level API uses the low-level function correctly
    date = Dates.DateTime(2000, 6, 21)
    lat = FT(45.0)

    # Using high-level API
    F_high, S_high, μ_high = daily_insolation(date, lat, param_set)

    # Using low-level functions manually
    Δt_years = Insolation.years_since_epoch(param_set, date)
    orb_params = Insolation.orbital_params(param_set)
    daily_θ, d =
        Insolation.daily_distance_zenith_angle(date, lat, orb_params, param_set)
    F_low, S_low, μ_low = insolation(daily_θ, d, param_set)

    # Should match exactly
    @test F_high == F_low
    @test S_high == S_low
    @test μ_high == μ_low
end

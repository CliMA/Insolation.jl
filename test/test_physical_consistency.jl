## Test physical consistency and bounds
rtol = 1e-3
atol = 1e-10

od = Insolation.OrbitalDataSplines()

@testset "Physical Bounds - Insolation Values" begin
    # Generate random dates and locations
    n_samples = 50

    for _ in 1:n_samples
        # Random date within reasonable range
        year = rand(1950:2050)
        month = rand(1:12)
        day = rand(1:28)  # Use 28 to avoid month-end issues
        hour = rand(0:23)
        date = Dates.DateTime(year, month, day, hour, 0, 0)

        # Random location
        lat = FT(rand(-90.0:0.1:90.0))
        lon = FT(rand(-180.0:0.1:180.0))

        # Test instantaneous insolation
        F, S, μ, ζ = insolation(date, lat, lon, param_set)

        # Insolation must be non-negative
        @test F >= 0

        # Solar flux must be positive
        @test S > 0

        # Cosine of zenith angle must be in [0, 1]
        @test 0 <= μ <= 1

        # Azimuth must be in [0, 2π]
        @test 0 <= ζ <= 2π

        # Insolation cannot exceed solar flux
        @test F <= S
    end
end

@testset "Physical Bounds - Daily Insolation" begin
    n_samples = 30

    for _ in 1:n_samples
        year = rand(1950:2050)
        month = rand(1:12)
        day = rand(1:28)
        date = Dates.DateTime(year, month, day)

        lat = FT(rand(-90.0:1.0:90.0))

        F, S, μ = daily_insolation(date, lat, param_set)

        # All values must be non-negative
        @test F >= 0
        @test S > 0

        # Daily averaged cosine must be in [0, 1]
        @test 0 <= μ <= 1

        # Daily insolation cannot exceed solar flux
        @test F <= S
    end
end

@testset "Physical Consistency - Solar Flux-Distance Relationship" begin
    # Solar flux should follow inverse square law with distance
    # S = S₀ * (d₀/d)²

    dates = [
        Dates.DateTime(2000, 1, 4),   # Near Earth's current perihelion
        Dates.DateTime(2000, 7, 4),   # Near Earth's current aphelion
    ]

    lat = FT(0.0)
    lon = FT(0.0)

    fluxes = Float64[]
    distances = Float64[]

    for date in dates
        Δt_years = Insolation.years_since_epoch(param_set, date)
        orb_params = Insolation.orbital_params(param_set)
        d, θ, ζ =
            Insolation.solar_geometry(date, lat, lon, orb_params, param_set)
        F, S, μ = insolation(θ, d, param_set)

        push!(fluxes, S)
        push!(distances, d)
    end

    # Test inverse square relationship
    # S₁/S₂ = (d₂/d₁)²
    ratio_flux = fluxes[1] / fluxes[2]
    ratio_dist_squared = (distances[2] / distances[1])^2

    @test ratio_flux ≈ ratio_dist_squared rtol = 1e-3
end

@testset "Physical Consistency - Zenith Angle at Noon" begin
    # At solar noon, zenith angle should be minimized
    date_base = Dates.DateTime(2000, 6, 21)
    lat = FT(45.0)
    lon = FT(0.0)

    # Sample throughout the day
    hours = 9:15
    zenith_angles = Float64[]

    for hour in hours
        date = date_base + Dates.Hour(hour)
        Δt_years = Insolation.years_since_epoch(param_set, date)
        orb_params = Insolation.orbital_params(param_set)
        d, θ, ζ =
            Insolation.solar_geometry(date, lat, lon, orb_params, param_set)
        push!(zenith_angles, θ)
    end

    # Minimum should occur around noon (index 4 corresponds to hour 12)
    min_idx = argmin(zenith_angles)
    @test min_idx in 3:5  # Allow for equation of time offset
end

@testset "Physical Consistency - Insolation = Flux × cos(SZA)" begin
    # Test the fundamental relationship: F = S × μ
    n_samples = 30

    for _ in 1:n_samples
        date = Dates.DateTime(
            rand(1950:2050),
            rand(1:12),
            rand(1:28),
            rand(0:23),
            0,
            0,
        )
        lat = FT(rand(-90.0:1.0:90.0))
        lon = FT(rand(-180.0:1.0:180.0))

        F, S, μ, ζ = insolation(date, lat, lon, param_set)

        # F should equal S × μ
        @test F ≈ S * μ rtol = 1e-5
    end
end

@testset "Physical Consistency - Earth-Sun Distance Bounds" begin
    # Earth-Sun distance should be within known bounds
    # Perihelion: ~147.1 million km
    # Aphelion: ~152.1 million km

    d0 = IP.orbit_semimaj(param_set)
    e = IP.eccentricity_epoch(param_set)

    # Perihelion distance
    d_peri = d0 * (1 - e)
    # Aphelion distance  
    d_aph = d0 * (1 + e)

    # Sample various dates
    dates = [Dates.DateTime(2000, m, 15) for m in 1:12]

    for date in dates
        Δt_years = Insolation.years_since_epoch(param_set, date)
        orb_params = Insolation.orbital_params(param_set)
        lat = FT(0.0)
        lon = FT(0.0)
        d, θ, ζ =
            Insolation.solar_geometry(date, lat, lon, orb_params, param_set)

        # Distance should be within orbital bounds
        @test d_peri <= d <= d_aph
    end
end

@testset "Physical Consistency - Time Continuity" begin
    # Insolation should change smoothly with time (no discontinuities)
    date_base = Dates.DateTime(2000, 6, 21, 12, 0, 0)
    lat = FT(45.0)
    lon = FT(0.0)

    # Sample at 1-minute intervals
    dt_minutes = 1
    n_samples = 60

    F_values = Float64[]
    for i in 0:(n_samples - 1)
        date = date_base + Dates.Minute(i * dt_minutes)
        F, _, _, _ = insolation(date, lat, lon, param_set)
        push!(F_values, F)
    end

    # Check that consecutive values don't differ too much
    for i in 1:(length(F_values) - 1)
        relative_change =
            abs(F_values[i + 1] - F_values[i]) / (F_values[i] + atol)
        @test relative_change < 0.01  # Less than 1% change per minute
    end
end

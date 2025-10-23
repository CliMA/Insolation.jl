atol = 1e-4
rtol = 1e-2
# relative tolerance for solar radiative flux/TSI comparison (which can differ by +/- 2 * eccentricity)
rtol_insol = 0.1

od = Insolation.OrbitalDataSplines()
## Test zero insolation at night
# sunrise at equator
date = Dates.DateTime(2020, 1, 1, 6, 0, 0)
lon, lat = [FT(0.0), FT(0.0)]
Δt_years = Insolation.years_since_epoch(param_set, date)
ϖ, γ, e = orbital_params(od, Δt_years)
orb_params = (FT(ϖ), FT(γ), FT(e))
d, θ, ζ = Insolation.solar_geometry(date, lat, lon, orb_params, param_set)
F, S, μ = insolation(θ, d, param_set)
@test F ≈ 0.0 atol = atol

# Test using high-level insolation API
F, S, μ, ζ = insolation(date, lat, lon, param_set, od; milankovitch = true)
@test S ≈ IP.tot_solar_irrad(param_set) rtol = rtol_insol
@test μ ≈ 0.0 atol = atol

# polar night NH 1
date = Dates.DateTime(2020, 12, 20, 11, 0, 0)
lon, lat = [FT(0.0), FT(80.0)]
Δt_years = Insolation.years_since_epoch(param_set, date)
ϖ, γ, e = orbital_params(od, Δt_years)
orb_params = (FT(ϖ), FT(γ), FT(e))
d, θ, ζ = Insolation.solar_geometry(date, lat, lon, orb_params, param_set)
F, S, μ = insolation(θ, d, param_set)
@test F ≈ 0.0 atol = atol

# Test using high-level insolation API
F, S, μ, ζ = insolation(date, lat, lon, param_set, od; milankovitch = true)
@test S ≈ IP.tot_solar_irrad(param_set) rtol = rtol_insol
@test μ ≈ 0.0 atol = atol

# polar night NH 2
date = Dates.DateTime(2020, 12, 20, 23, 0, 0)
Δt_years = Insolation.years_since_epoch(param_set, date)
ϖ, γ, e = orbital_params(od, Δt_years)
orb_params = (FT(ϖ), FT(γ), FT(e))
d, θ, ζ = Insolation.solar_geometry(date, lat, lon, orb_params, param_set)
F, S, μ = insolation(θ, d, param_set)
@test F ≈ 0.0 atol = atol

# Test using high-level insolation API
F, S, μ, ζ = insolation(date, lat, lon, param_set, od; milankovitch = true)
@test S ≈ IP.tot_solar_irrad(param_set) rtol = rtol_insol
@test μ ≈ 0.0 atol = atol

## Test symmetry of insolation at equinox
nlats = 181
date = Dates.DateTime(2021, 3, 20, 9, 37, 0) # vernal equinox 2021
l_arr = FT.(collect(range(-90, stop = 90, length = nlats)))
F_arr = zeros(nlats)
for (i, lat) in enumerate(l_arr)
    local Δt_years = Insolation.years_since_epoch(param_set, date)
    local ϖ, γ, e = orbital_params(od, Δt_years)
    local orb_params = (FT(ϖ), FT(γ), FT(e))
    local daily_θ, d =
        Insolation.daily_distance_zenith_angle(date, lat, orb_params, param_set)
    local F, S, μ = insolation(daily_θ, d, param_set)
    F_arr[i] = F
end
F_NH = sort(F_arr[l_arr .>= 0])
F_SH = sort(F_arr[l_arr .<= 0])
@test F_NH ≈ F_SH rtol = rtol

## Test globally averaged insolation ≈ TSI
ndays, nlats = [365, 361]
d_arr = Array{Int}(round.(collect(range(0, stop = 365, length = ndays))))
l_arr = FT.(collect(range(-90, stop = 90, length = nlats)))
F_arr = zeros(ndays, nlats)

for (i, d) in enumerate(d_arr)
    for (j, lat) in enumerate(l_arr)
        local datei = Dates.DateTime(2000, 1, 1) + Dates.Day(d)
        local Δt_years = Insolation.years_since_epoch(param_set, datei)
        local ϖ, γ, e = orbital_params(od, Δt_years)
        local orb_params = (FT(ϖ), FT(γ), FT(e))
        local daily_θ, dist = Insolation.daily_distance_zenith_angle(
            datei,
            lat,
            orb_params,
            param_set,
        )
        local F, S, μ = insolation(daily_θ, dist, param_set)
        F_arr[i, j] = F
    end
end

zonal_mean_insol = mean(F_arr, dims = 1)
area_fac = abs.(cosd.(l_arr))
global_mean_insol = sum(zonal_mean_insol * area_fac) / sum(area_fac)
@test global_mean_insol ≈ IP.tot_solar_irrad(param_set) / 4 rtol = rtol

## Test invariance of zonal-mean insolation under rotation of ϖ
ϖ0 = IP.lon_perihelion_epoch(param_set)
param_set = IP.InsolationParameters(FT, (; lon_perihelion_epoch = ϖ0 + π))

for (i, d) in enumerate(d_arr)
    for (j, lat) in enumerate(l_arr)
        local datei = Dates.DateTime(2000, 1, 1) + Dates.Day(d)
        local Δt_years = Insolation.years_since_epoch(param_set, datei)
        local ϖ, γ, e = orbital_params(od, Δt_years)
        local orb_params = (FT(ϖ), FT(γ), FT(e))
        local daily_θ, dist = Insolation.daily_distance_zenith_angle(
            datei,
            lat,
            orb_params,
            param_set,
        )
        local F, S, μ = insolation(daily_θ, dist, param_set)
        F_arr[i, j] = F
    end
end

zonal_mean_insol_rotate = mean(F_arr, dims = 1)
@test zonal_mean_insol_rotate ≈ zonal_mean_insol rtol = rtol

param_set = IP.InsolationParameters(FT, (; lon_perihelion_epoch = ϖ0))

## Test simplified insolation method (without OrbitalDataSplines)
# Test that simplified method works for various dates and locations
date = Dates.DateTime(2000, 6, 21, 12, 0, 0)  # Summer solstice
lon, lat = [FT(0.0), FT(45.0)]  # 45°N at Greenwich meridian

# Test simplified method (no OrbitalDataSplines needed)
F_simple, S_simple, μ_simple, ζ_simple = insolation(date, lat, lon, param_set)
@test S_simple > 0.0
@test μ_simple > 0.0
@test μ_simple <= 1.0

# Test full method with OrbitalDataSplines
F_full, S_full, μ_full, ζ_full =
    insolation(date, lat, lon, param_set, od; milankovitch = true)

# For year 2000 (near epoch), both methods should give very similar results
@test S_simple ≈ S_full rtol = 0.01
@test μ_simple ≈ μ_full rtol = 0.01

# Test at night - both should give μ = 0
date_night = Dates.DateTime(2000, 6, 21, 0, 0, 0)  # Midnight
F_night, S_night, μ_night, ζ_night = insolation(date_night, lat, lon, param_set)
@test μ_night ≈ 0.0 atol = atol

# Test at equator at noon near equinox - should have high μ
date_equinox = Dates.DateTime(2000, 3, 20, 12, 0, 0)
lon_eq, lat_eq = [FT(0.0), FT(0.0)]
F_eq, S_eq, μ_eq, ζ_eq = insolation(date_equinox, lat_eq, lon_eq, param_set)
@test μ_eq > 0.95  # Sun nearly overhead at equator during equinox noon

# Test error when milankovitch=true but orbital_data is nothing
@test_throws ErrorException insolation(
    date,
    lat,
    lon,
    param_set,
    nothing;
    milankovitch = true,
)

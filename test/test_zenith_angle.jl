rtol = 1e-2
od = Insolation.OrbitalDataSplines()

# sunrise at equator
date = Dates.DateTime(2020, 2, 20, 6, 11, 0)
lon, lat = [FT(0.0), FT(0.0)]
Δt_years = Insolation.years_since_epoch(param_set, date)
ϖ, γ, e = orbital_params(od, Δt_years)
orb_params = (FT(ϖ), FT(γ), FT(e))
(; d, θ, ζ) = solar_geometry(date, lat, lon, orb_params, param_set)
@test θ ≈ π / 2 rtol = rtol

# solar noon at equator
date = Dates.DateTime(2020, 2, 20, 12, 14, 0)
Δt_years = Insolation.years_since_epoch(param_set, date)
ϖ, γ, e = orbital_params(od, Δt_years)
orb_params = (FT(ϖ), FT(γ), FT(e))
(; d, θ, ζ) = solar_geometry(date, lat, lon, orb_params, param_set)
@test ζ ≈ 3π / 2 rtol = rtol

# sunset at equator
date = Dates.DateTime(2020, 2, 20, 18, 17, 0)
Δt_years = Insolation.years_since_epoch(param_set, date)
ϖ, γ, e = orbital_params(od, Δt_years)
orb_params = (FT(ϖ), FT(γ), FT(e))
(; d, θ, ζ) = solar_geometry(date, lat, lon, orb_params, param_set)
@test θ ≈ π / 2 rtol = rtol

# solar noon at equator, eot correction = false
date = Dates.DateTime(2020, 2, 20, 12, 0, 0)
Δt_years = Insolation.years_since_epoch(param_set, date)
ϖ, γ, e = orbital_params(od, Δt_years)
orb_params = (FT(ϖ), FT(γ), FT(e))
(; d, θ, ζ) = solar_geometry(date, lat, lon, orb_params, param_set, eot_correction = false)
@test ζ ≈ 3π / 2 rtol = rtol

# sunset at equator, eot correction = false
date = Dates.DateTime(2020, 2, 20, 18, 0, 0)
Δt_years = Insolation.years_since_epoch(param_set, date)
ϖ, γ, e = orbital_params(od, Δt_years)
orb_params = (FT(ϖ), FT(γ), FT(e))
(; d, θ, ζ) = solar_geometry(date, lat, lon, orb_params, param_set, eot_correction = false)
@test θ ≈ π / 2 rtol = rtol

## Test Polar Night
# polar night NH 1
date = Dates.DateTime(2020, 12, 20, 11, 0, 0)
lon, lat = [FT(0.0), FT(80.0)]
Δt_years = Insolation.years_since_epoch(param_set, date)
ϖ, γ, e = orbital_params(od, Δt_years)
orb_params = (FT(ϖ), FT(γ), FT(e))
(; d, θ, ζ) = solar_geometry(date, lat, lon, orb_params, param_set)
@test θ > π / 2

# polar night NH 2
date = Dates.DateTime(2020, 12, 20, 23, 0, 0)
Δt_years = Insolation.years_since_epoch(param_set, date)
ϖ, γ, e = orbital_params(od, Δt_years)
orb_params = (FT(ϖ), FT(γ), FT(e))
(; d, θ, ζ) = solar_geometry(date, lat, lon, orb_params, param_set)
@test θ > π / 2

# polar night SH 1
date = Dates.DateTime(2020, 6, 20, 11, 0, 0)
lon, lat = [FT(0.0), FT(-80.0)]
Δt_years = Insolation.years_since_epoch(param_set, date)
ϖ, γ, e = orbital_params(od, Δt_years)
orb_params = (FT(ϖ), FT(γ), FT(e))
(; d, θ, ζ) = solar_geometry(date, lat, lon, orb_params, param_set)
@test θ > π / 2

# polar night SH 2
date = Dates.DateTime(2020, 6, 20, 23, 0, 0)
Δt_years = Insolation.years_since_epoch(param_set, date)
ϖ, γ, e = orbital_params(od, Δt_years)
orb_params = (FT(ϖ), FT(γ), FT(e))
(; d, θ, ζ) = solar_geometry(date, lat, lon, orb_params, param_set)
@test θ > π / 2

## Test Distance
date = Dates.DateTime(2000, 3, 22, 0, 0, 0)
Δt_years = Insolation.years_since_epoch(param_set, date)
ϖ, γ, e = orbital_params(od, Δt_years)
orb_params = (FT(ϖ), FT(γ), FT(e))
(; d, θ, ζ) = solar_geometry(date, lat, lon, orb_params, param_set)
@test d ≈ IP.orbit_semimaj(param_set) rtol = rtol

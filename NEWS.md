# Insolation.jl Release Notes

## v1.2.0

### Breaking changes

- Removed the `orbit_semimaj` field from `InsolationParameters` (and the
  `length_orbit_semi_major` ClimaParams mapping). The total solar irradiance
  `tot_solar_irrad` is defined at the *mean* orbital distance (the semi-major
  axis), so the absolute orbit size cancels out of the insolation and is no longer
  needed. As a result, the planet-star distance `d` returned by `solar_geometry`,
  `daily_distance_zenith_angle`, and `planet_star_distance` (and accepted by the
  low-level `solar_flux`/`insolation(θ, d, …)` methods) is now **dimensionless**,
  expressed in units of the semi-major axis (≈ AU for Earth) rather than meters.
  Multiply by the semi-major axis to recover a physical distance. The flux is
  `S0 / d^2`. This does not affect `F`, `S`, `μ`, or `ζ`.

### Bug fixes

- The equation of time is now computed from the exact equatorial projection of the
  Sun's ecliptic longitude (`Δη = L − α`, with `α = atan2(cosγ sinLₛ, cosLₛ)`)
  instead of the small-obliquity expansion `tan²(γ/2)`. It reduces to the previous
  formula for small tilt but remains valid for any obliquity (e.g., Uranus-like).
- Fixed an interval-selection bug in `TSIDataSpline`'s `evaluate` that produced a
  small non-physical discontinuity in the interpolated total solar irradiance for
  dates after the 15th of a month but before noon.
- The Milankovitch (paleoclimate) path now indexes the Laskar et al. (2004) orbital
  tables in Julian years, matching the tables' native time axis, instead of
  anomalistic years. This slightly changes time-varying orbital-parameter results,
  most noticeably at deep-time extremes.
- `solar_geometry` and `daily_distance_zenith_angle` now convert all inputs to
  `eltype(param_set)` internally, removing silent floating-point type mixing (e.g.,
  Float32 parameter sets no longer leak Float64 into the result).

### Enhancements

- Added `julian_years_since_epoch` for computing time since epoch in Julian years
  (the unit of the Laskar tables), complementing the anomalistic-year
  `years_since_epoch`.
- Made the solar azimuth calculation robust at high declinations by removing a
  `tan(δ)` singularity (relevant for extreme-obliquity configurations).
- Documentation audit: corrected physical descriptions and stale examples,
  standardized docstrings, and added a type-level docstring for
  `InsolationParameters`. Documented that the `true_anomaly` series (and the
  quantities derived from it) is a low-eccentricity approximation that should be
  replaced by an exact Kepler solver for highly eccentric orbits.
- Clarified that `tot_solar_irrad` is the total solar irradiance at the *mean*
  orbital distance (the semi-major axis), consistent with ClimaParams (earlier docs
  incorrectly described this as "at 1 au").
- Documented that the `TSIDataSpline` is the Sun's irradiance at 1 au and is
  therefore meaningful for Earth only: in the flux calculation it is taken as the
  irradiance at the planet's mean orbital distance, which holds only because
  Earth's semi-major axis is ≈ 1 au (the package does not rescale it for other
  bodies).

### Notes

- The signatures of `solar_geometry` and `daily_distance_zenith_angle` were widened
  to accept any tuple of real-valued orbital parameters; existing calls are
  unaffected.
- Several bug fixes change numerical results slightly for deep-time Milankovitch
  calculations and for time-varying TSI (1e-3 W/m2 interpolation errors corrected)

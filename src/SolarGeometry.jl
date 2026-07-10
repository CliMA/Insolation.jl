export solar_geometry, daily_distance_zenith_angle

"""
    mean_anomaly(Δt_years::FT, param_set::IP.AIP) where {FT}

Calculate the mean anomaly at a given time since epoch.

The mean anomaly is the angle the planet would have traveled from perihelion
if it moved in a circular orbit at constant angular velocity.

# Arguments
- `Δt_years::FT`: Time since epoch [anomalistic years]
- `param_set::IP.AIP`: Parameter struct containing `mean_anom_epoch`

# Returns
- `MA`: Mean anomaly [radians]
"""
function mean_anomaly(Δt_years::FT, param_set::IP.AIP) where {FT}
    M0 = IP.mean_anom_epoch(param_set)
    MA = mod(FT(2π) * Δt_years + M0, FT(2π))
    return MA
end

"""
    true_anomaly(MA::FT, e::FT) where {FT <: Real}

Calculate the true anomaly from the mean anomaly.

The true anomaly is the actual angular distance of the planet from perihelion
along its orbital path. This function uses a series expansion (the "equation of
the center") that is accurate to third order in the eccentricity `e`, with an
error of O(e⁴) (see Fitzpatrick (2012), Appendix A.10).

!!! warning "Low-eccentricity approximation"
    The series is intended for nearly circular orbits such as Earth's
    (`e ≈ 0.0167`). Its error grows rapidly with eccentricity (exceeding a few
    degrees near `e ≈ 0.5`), so for highly eccentric orbits it should be replaced
    by an exact solution of Kepler's equation (e.g., a Newton–Raphson iteration
    for the eccentric anomaly). This approximation affects every quantity derived
    from the true anomaly (declination, distance, and the equation of time).

# Arguments
- `MA::FT`: Mean anomaly [radians]
- `e::FT`: Orbital eccentricity [unitless]

# Returns
- `TA`: True anomaly [radians]
"""
function true_anomaly(MA::FT, e::FT) where {FT<:Real}
    # Series expansion for true anomaly
    TA =
        MA +
        (2 * e - FT(1 / 4) * e^3) * sin(MA) +
        FT(5 / 4) * e^2 * sin(2 * MA) +
        FT(13 / 12) * e^3 * sin(3 * MA)
    return mod(TA, FT(2π))
end

"""
    solar_longitude(TA::FT, ϖ::FT) where {FT <: Real}

Calculate the solar longitude (ecliptic longitude of the Sun).

The solar longitude is the angular distance of the planet along its orbital
path, measured from the vernal equinox. It is the sum of the true anomaly
(angle from perihelion) and the longitude of perihelion.

# Arguments
- `TA::FT`: True anomaly [radians]
- `ϖ::FT`: Longitude of perihelion [radians]

# Returns
- `SL`: Solar longitude [radians]
"""
function solar_longitude(TA::FT, ϖ::FT) where {FT<:Real}
    SL = mod(TA + ϖ, FT(2π))
    return SL
end

"""
    hour_angle(
        date::DateTime,
        λ::FT,
        MA::FT,
        (ϖ, γ, e)::Tuple{FT, FT, FT};
        eot_correction::Bool = true,
    ) where {FT}

Calculate the hour angle at a given time and longitude.

The hour angle is zero at local solar noon and increases with time.
Optionally applies the equation of time correction to account for 
the difference between apparent and mean solar time.

# Arguments
- `date::DateTime`: Current date and time
- `λ::FT`: Longitude [radians]
- `MA::FT`: Mean anomaly [radians]
- `(ϖ, γ, e)::Tuple{FT, FT, FT}`: Orbital parameters tuple containing:
  - `ϖ`: Longitude of perihelion [radians]
  - `γ`: Obliquity (axial tilt) [radians]
  - `e`: Orbital eccentricity [unitless]
- `eot_correction::Bool`: (default true) Apply equation of time correction

# Returns
- `η`: Hour angle [radians]
"""
function hour_angle(
    date::DateTime,
    λ::FT,
    MA::FT,
    (ϖ, γ, e)::Tuple{FT,FT,FT};
    eot_correction::Bool = true,
) where {FT}
    # Equation of time correction 
    Δη_eot = equation_of_time(MA, (ϖ, γ, e))
    Δη = ifelse(eot_correction, Δη_eot, FT(0))

    time_of_day = FT(mod(datetime2julian(date), 1))  # fractional day [0, 1)
    η_prime_uncorrected = FT(2π) * time_of_day       # uncorrected hour angle [radians]
    η_prime = mod(η_prime_uncorrected + Δη, FT(2π))

    # Hour angle at given longitude [given in radians]
    η = mod(η_prime + λ, FT(2π))
    return η
end

"""
    equation_of_time(MA::FT, (ϖ, γ, e)::Tuple{FT, FT, FT}) where {FT <: Real}

Calculate the equation of time correction for the hour angle.

The equation of time accounts for the difference between apparent solar time
(based on the actual Sun's position in the sky) and mean solar time (based on
a fictitious mean Sun moving at constant speed). This difference arises from
two effects:
1. The elliptical orbit (eccentricity e)
2. The axial tilt (obliquity γ)

It is computed as the difference between the mean longitude ``L = M + \\varpi``
and the right ascension ``\\alpha`` of the true Sun. The right ascension is
obtained from the exact projection of the ecliptic longitude onto the equatorial
plane, ``\\alpha = \\mathrm{atan2}(\\cos\\gamma \\sin L_s, \\cos L_s)`` with
``L_s`` the solar longitude, so the obliquity contribution is exact for any
``\\gamma`` (it does not rely on a small-tilt expansion). The eccentricity
contribution enters through the true anomaly and is therefore accurate to the
same order in `e` as [`true_anomaly`](@ref).

# Arguments
- `MA::FT`: Mean anomaly [radians]
- `(ϖ, γ, e)::Tuple{FT, FT, FT}`: Orbital parameters tuple containing:
  - `ϖ`: Longitude of perihelion [radians]
  - `γ`: Obliquity (axial tilt) [radians]
  - `e`: Orbital eccentricity [unitless]

# Returns
- `Δη`: Hour angle correction [radians]
"""
function equation_of_time(MA::FT, (ϖ, γ, e)::Tuple{FT,FT,FT}) where {FT<:Real}
    # Mean longitude of the fictitious mean Sun and true ecliptic (solar) longitude
    L = MA + ϖ
    L_s = true_anomaly(MA, e) + ϖ
    # Right ascension of the true Sun from the exact equatorial projection (valid
    # for any obliquity γ, including γ > π/2)
    α = atan(cos(γ) * sin(L_s), cos(L_s))
    # Equation of time as an hour-angle correction: mean longitude − right ascension
    Δη = L - α
    return mod(Δη + FT(π), FT(2π)) - FT(π)
end

"""
    planet_star_distance(TA::FT, e::FT) where {FT <: Real}

Calculate the planet-star distance at a given true anomaly, in units of the
semi-major axis (the mean orbital distance).

The distance varies due to the planet's elliptical orbit, being shortest at
perihelion and longest at aphelion. The calculation uses the orbit equation
for an ellipse. The result is normalized by the semi-major axis ``a``, so it is
the dimensionless ratio ``d/a = (1-e^2)/(1+e\\cos A)``, which ranges from `1-e`
(perihelion) to `1+e` (aphelion). The absolute orbit size is irrelevant to the
insolation, so it is not carried; multiply by the semi-major axis to recover a
physical distance.

# Arguments
- `TA::FT`: True anomaly [radians]
- `e::FT`: Orbital eccentricity [unitless]

# Returns
- `d`: Planet-star distance, in units of the semi-major axis [unitless]
"""
function planet_star_distance(TA::FT, e::FT) where {FT<:Real}
    d = (1 - e^2) / (1 + e * cos(TA))
    return d
end

"""
    years_since_epoch(
        param_set::IP.InsolationParameters{FT},
        date::DateTime,
    ) where {FT}

Calculate the time elapsed since epoch (typically J2000) in anomalistic
years (the time from perihelion to perihelion).

Converts the time difference between two dates from Julian days to
anomalistic years, which is the natural time unit for orbital calculations.

# Arguments
- `param_set::IP.InsolationParameters{FT}`: Parameter struct
- `date::DateTime`: Current date and time

# Returns
- `Δt_years`: Time since epoch [anomalistic years]

See also [`julian_years_since_epoch`](@ref), which uses Julian years (the time unit
of the Laskar et al. (2004) orbital tables).
"""
function years_since_epoch(
    param_set::IP.InsolationParameters{FT},
    date::DateTime,
) where {FT}
    (; epoch, year_anom, day) = param_set
    days_per_year = year_anom / day
    return FT(datetime2julian(date) - datetime2julian(epoch)) / days_per_year
end

"""
    julian_years_since_epoch(
        param_set::IP.InsolationParameters{FT},
        date::DateTime,
    ) where {FT}

Calculate the time elapsed since epoch (typically J2000) in Julian years
(a Julian year is exactly 365.25 days of 86400 s each).

This is the time unit of the Laskar et al. (2004) orbital-parameter tables, so it
is used to index [`OrbitalDataSplines`](@ref). It differs slightly from
[`years_since_epoch`](@ref), which returns *anomalistic* years (the natural unit
for the mean anomaly); the two differ by the ratio of the anomalistic year to the
Julian year (~0.0026% for Earth).

# Arguments
- `param_set::IP.InsolationParameters{FT}`: Parameter struct
- `date::DateTime`: Current date and time

# Returns
- `Δt_years`: Time since epoch [Julian years]
"""
function julian_years_since_epoch(
    param_set::IP.InsolationParameters{FT},
    date::DateTime,
) where {FT}
    epoch = IP.epoch(param_set)
    days_per_julian_year = FT(365.25)
    return FT(datetime2julian(date) - datetime2julian(epoch)) / days_per_julian_year
end

"""
    distance_declination_mean_anomaly(
        Δt_years::FT,
        (ϖ, γ, e)::Tuple{FT, FT, FT},
        param_set::IP.AIP,
    ) where {FT}

Compute the planet-star distance, solar declination angle, and mean anomaly.

This function calculates key astronomical parameters from orbital elements.
The declination determines the subsolar latitude, while the planet-star distance
varies due to orbital eccentricity. The mean anomaly is returned for use in
hour angle calculations.

# Arguments
- `Δt_years::FT`: Time since epoch [anomalistic years]
- `(ϖ, γ, e)::Tuple{FT, FT, FT}`: Orbital parameters tuple containing:
  - `ϖ`: Longitude of perihelion [radians]
  - `γ`: Obliquity (axial tilt) [radians]
  - `e`: Orbital eccentricity [unitless]
- `param_set::IP.AIP`: Parameter struct

# Returns
- `d`: Planet-star distance, in units of the semi-major axis [unitless]
- `δ`: Solar declination angle [radians]
- `MA`: Mean anomaly [radians]
"""
function distance_declination_mean_anomaly(
    Δt_years::FT,
    (ϖ, γ, e)::Tuple{FT,FT,FT},
    param_set::IP.AIP,
) where {FT}
    # Mean anomaly at current time
    MA = mean_anomaly(Δt_years, param_set)

    # True anomaly [radians]
    TA = true_anomaly(MA, e)

    # Solar longitude [radians]
    SL = solar_longitude(TA, ϖ)

    # Declination angle [radians], in [-π/2, π/2]
    δ = asin(sin(γ) * sin(SL))

    # Planet-star distance, in units of the semi-major axis
    d = planet_star_distance(TA, e)

    return d, δ, MA
end

"""
    solar_geometry(
        date::DateTime,
        latitude::Real,
        longitude::Real,
        orb_params::Tuple{<:Real, <:Real, <:Real},
        param_set::AIP;
        eot_correction::Bool = true,
    )

Calculate the planet-star distance, solar zenith angle, and azimuth angle.

This is a high-level function that combines all necessary astronomical
calculations to determine the position of the star in the sky (zenith and
azimuth angles) and the planet-star distance at a specific time and location.

All real-valued inputs are converted to `eltype(param_set)` internally, so the
computation is carried out consistently in the parameter set's floating-point type.

# Arguments
- `date::DateTime`: Current date and time
- `latitude::Real`: Latitude [degrees]
- `longitude::Real`: Longitude [degrees]
- `orb_params::Tuple{<:Real, <:Real, <:Real}`: Orbital parameters tuple `(ϖ, γ, e)`:
  - `ϖ`: Longitude of perihelion [radians]
  - `γ`: Obliquity (axial tilt) [radians]
  - `e`: Orbital eccentricity [unitless]
- `param_set::AIP`: Parameter struct
- `eot_correction::Bool`: (default true) Apply equation of time correction

# Returns
A `NamedTuple` with fields:
- `d`: Planet-Sun distance, in units of the semi-major axis [unitless]
- `θ`: Solar zenith angle [radians]
- `ζ`: Solar azimuth angle [radians], 0 = due East, increasing CCW
"""
function solar_geometry(
    date::DateTime,
    latitude::Real,
    longitude::Real,
    orb_params::Tuple{<:Real,<:Real,<:Real},
    param_set::AIP;
    eot_correction = true,
)
    FT = eltype(param_set)
    # Convert all inputs to the parameter set's float type for type consistency
    ϖ, γ, e = FT.(orb_params)
    ϕ = FT(deg2rad(latitude))
    λ = FT(deg2rad(longitude))

    # Get time since epoch in anomalistic years
    Δt_years = years_since_epoch(param_set, date)

    # Get distance, declination, and mean anomaly
    d, δ, MA = distance_declination_mean_anomaly(Δt_years, (ϖ, γ, e), param_set)

    # Get hour angle at given longitude
    η = hour_angle(date, λ, MA, (ϖ, γ, e); eot_correction)

    # Solar zenith angle [radians], in [0, π] (argument clamped for acos safety)
    θ = acos(clamp(cos(ϕ) * cos(δ) * cos(η) + sin(ϕ) * sin(δ), FT(-1), FT(1)))

    # Solar azimuth angle: ζ = 0 when due E and increasing CCW
    # ζ = 3π/2 (due S) when η=0 at local solar noon.
    # The arguments are multiplied through by cos(δ) ≥ 0 to avoid the tan(δ)
    # singularity at high declinations (e.g., extreme obliquity).
    ζ = mod(
        FT(3π / 2) - atan(cos(δ) * sin(η), cos(δ) * cos(η) * sin(ϕ) - sin(δ) * cos(ϕ)),
        FT(2π),
    )

    return (; d, θ, ζ)
end

"""
    daily_distance_zenith_angle(
        date::DateTime,
        latitude::Real,
        orb_params::Tuple{<:Real, <:Real, <:Real},
        param_set::IP.AIP,
    )

Calculate the effective zenith angle for diurnally averaged insolation and the
planet-star distance.

Return the effective zenith angle corresponding to the diurnally averaged
insolation (averaging the cosine of the zenith angle over 24 hours) and the
planet-star distance for a given date and latitude.

All real-valued inputs are converted to `eltype(param_set)` internally, so the
computation is carried out consistently in the parameter set's floating-point type.

# Arguments
- `date::DateTime`: Current date
- `latitude::Real`: Latitude [degrees]
- `orb_params::Tuple{<:Real, <:Real, <:Real}`: Orbital parameters tuple `(ϖ, γ, e)`:
  - `ϖ`: Longitude of perihelion [radians]
  - `γ`: Obliquity (axial tilt) [radians]
  - `e`: Orbital eccentricity [unitless]
- `param_set::IP.AIP`: Parameter struct

# Returns
A `NamedTuple` with fields:
- `daily_θ`: Effective solar zenith angle [radians]
- `d`: Planet-star distance, in units of the semi-major axis [unitless]
"""
function daily_distance_zenith_angle(
    date::DateTime,
    latitude::Real,
    orb_params::Tuple{<:Real,<:Real,<:Real},
    param_set::IP.AIP,
)
    FT = eltype(param_set)
    # Convert all inputs to the parameter set's float type for type consistency
    ϖ, γ, e = FT.(orb_params)
    ϕ = FT(deg2rad(latitude))

    Δt_years = years_since_epoch(param_set, date)

    # Get distance and declination
    d, δ, _ = distance_declination_mean_anomaly(Δt_years, (ϖ, γ, e), param_set)

    # Sunrise/sunset hour angle. Clamping the argument to [-1, 1] also encodes the
    # polar limits in a single branchless expression: clamp gives acos(-1) = π for
    # polar day (tanϕ·tanδ ≥ 1) and acos(1) = 0 for polar night (tanϕ·tanδ ≤ -1).
    ηd = acos(clamp(-tan(ϕ) * tan(δ), FT(-1), FT(1)))

    # Effective zenith angle to get diurnally averaged insolation
    # (i.e., averaging the cosine of the zenith angle). The argument is clamped
    # to [-1, 1] to guard against floating-point overshoot in acos.
    daily_cosθ =
        clamp(FT(1 / π) * (ηd * sin(ϕ) * sin(δ) + cos(ϕ) * cos(δ) * sin(ηd)), FT(-1), FT(1))
    daily_θ = acos(daily_cosθ)

    return (; daily_θ, d)
end

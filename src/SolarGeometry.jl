export solar_geometry, daily_distance_zenith_angle

"""
    mean_anomaly(Δt_years::FT, param_set::IP.AIP) where {FT}

Calculates the mean anomaly at a given time since epoch.

The mean anomaly is the angle the planet would have traveled from perihelion
if it moved in a circular orbit at constant angular velocity.

# Arguments
- `Δt_years::FT`: Time since epoch [years]
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

Calculates the true anomaly from the mean anomaly.

The true anomaly is the actual angular distance of the planet from perihelion
along its orbital path. This function uses a series expansion accurate to
O(e⁴) where e is the eccentricity (see Fitzpatrick (2012), Appendix A.10).

# Arguments
- `MA::FT`: Mean anomaly [radians]
- `e::FT`: Orbital eccentricity [unitless]

# Returns
- `TA`: True anomaly [radians]
"""
function true_anomaly(MA::FT, e::FT) where {FT <: Real}
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

Calculates the solar longitude (ecliptic longitude of the Sun).

The solar longitude is the angular distance of the planet along its orbital 
path, measured from vernal equinox. It is the sum of the true anomaly 
(angle from perihelion) and the longitude of perihelion.

# Arguments
- `TA::FT`: True anomaly [radians]
- `ϖ::FT`: Longitude of perihelion [radians]

# Returns
- `SL`: Solar longitude [radians]
"""
function solar_longitude(TA::FT, ϖ::FT) where {FT <: Real}
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

Calculates the hour angle at a given time and longitude.

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
    (ϖ, γ, e)::Tuple{FT, FT, FT};
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

Calculates the equation of time correction for the hour angle.

The equation of time accounts for the difference between apparent solar time 
(based on the actual Sun's position in the sky) and mean solar time (based on 
a fictitious mean Sun moving at constant speed). This difference arises from 
two effects:
1. Earth's elliptical orbit (eccentricity e)
2. Earth's axial tilt (obliquity γ)

# Arguments
- `MA::FT`: Mean anomaly [radians]
- `(ϖ, γ, e)::Tuple{FT, FT, FT}`: Orbital parameters tuple containing:
  - `ϖ`: Longitude of perihelion [radians]
  - `γ`: Obliquity (axial tilt) [radians]
  - `e`: Orbital eccentricity [unitless]

# Returns
- `Δη`: Hour angle correction [radians]
"""
function equation_of_time(
    MA::FT,
    (ϖ, γ, e)::Tuple{FT, FT, FT},
) where {FT <: Real}
    Δη = -2 * e * sin(MA) + tan(γ / 2)^2 * sin(2 * (MA + ϖ))
    return mod(Δη + FT(π), FT(2π)) - FT(π)
end

"""
    planet_star_distance(TA::FT, e::FT, param_set::IP.AIP) where {FT <: Real}

Calculates the distance between planet (Earth) and star (Sun) at a given 
true anomaly.

The distance varies due to the planet's elliptical orbit, being shortest at
perihelion and longest at aphelion. The calculation uses the orbit equation
for an ellipse.

# Arguments
- `TA::FT`: True anomaly [radians]
- `e::FT`: Orbital eccentricity [unitless]
- `param_set::IP.AIP`: Parameter struct containing `orbit_semimaj`

# Returns
- `d`: Planet-star distance [m]
"""
function planet_star_distance(
    TA::FT,
    e::FT,
    param_set::IP.AIP,
) where {FT <: Real}
    d0 = IP.orbit_semimaj(param_set)
    d = d0 * (1 - e^2) / (1 + e * cos(TA))
    return d
end

"""
    years_since_epoch(
        param_set::IP.InsolationParameters{FT},
        date::DateTime,
    ) where {FT}

Calculates the time elapsed since epoch (typically J2000) in anomalistic 
years (the time from perihelion to perihelion).

Converts the time difference between two dates from Julian days to
anomalistic years, which is the natural time unit for orbital calculations.

# Arguments
- `param_set::IP.InsolationParameters{FT}`: Parameter struct
- `date::DateTime`: Current date and time

# Returns
- `Δt_years`: Time since epoch [anomalistic years]
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
    distance_declination_mean_anomaly(
        Δt_years::FT,
        (ϖ, γ, e)::Tuple{FT, FT, FT},
        param_set::IP.AIP,
    ) where {FT}

Computes planet-star distance, solar declination angle, and mean anomaly.

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
- `d`: Planet-star distance [m]
- `δ`: Solar declination angle [radians]
- `MA`: Mean anomaly [radians]
"""
function distance_declination_mean_anomaly(
    Δt_years::FT,
    (ϖ, γ, e)::Tuple{FT, FT, FT},
    param_set::IP.AIP,
) where {FT}
    # Mean anomaly at current time
    MA = mean_anomaly(Δt_years, param_set)

    # True anomaly [radians]
    TA = true_anomaly(MA, e)

    # Solar longitude [radians]
    SL = solar_longitude(TA, ϖ)

    # Declination angle [radians]
    δ = mod(asin(sin(γ) * sin(SL)), FT(2π))

    # Planet-star distance [m] 
    d = planet_star_distance(TA, e, param_set)

    return d, δ, MA
end

"""
    solar_geometry(
        date::DateTime,
        latitude::Real,
        longitude::Real,
        (ϖ, γ, e)::Tuple{FT, FT, FT},
        param_set::AIP;
        eot_correction::Bool = true,
    ) where {FT}

Calculates planet-star distance, solar zenith angle, and azimuthal angle.

This is a high-level function that combines all necessary astronomical
calculations to determine the planet's position and distance from the star
at a specific time and location.

# Arguments
- `date::DateTime`: Current date and time
- `latitude::Real`: Latitude [degrees]
- `longitude::Real`: Longitude [degrees]
- `(ϖ, γ, e)::Tuple{FT, FT, FT}`: Orbital parameters tuple containing:
  - `ϖ`: Longitude of perihelion [radians]
  - `γ`: Obliquity (axial tilt) [radians]
  - `e`: Orbital eccentricity [unitless]
- `param_set::AIP`: Parameter struct
- `eot_correction::Bool`: (default true) Apply equation of time correction

# Returns
- `d`: Planet-Sun distance [m]
- `θ`: Solar zenith angle [radians]
- `ζ`: Solar azimuth angle [radians], 0 = due East, increasing CCW
"""
function solar_geometry(
    date::DateTime,
    latitude::Real,
    longitude::Real,
    (ϖ, γ, e)::Tuple{FT, FT, FT},
    param_set::AIP;
    eot_correction = true,
) where {FT}
    ϕ = eltype(param_set)(deg2rad(latitude))
    λ = eltype(param_set)(deg2rad(longitude))

    # Get time since epoch in anomalistic years
    Δt_years = years_since_epoch(param_set, date)

    # Get distance, declination, and mean anomaly
    d, δ, MA = distance_declination_mean_anomaly(Δt_years, (ϖ, γ, e), param_set)

    # Get hour angle at given longitude
    η = hour_angle(date, λ, MA, (ϖ, γ, e); eot_correction)

    # Solar zenith angle [radians]
    θ = mod(
        acos(
            max(FT(-1), min(FT(1), cos(ϕ) * cos(δ) * cos(η) + sin(ϕ) * sin(δ))),
        ),
        FT(2π),
    )

    # Solar azimuth angle: ζ = 0 when due E and increasing CCW
    # ζ = 3π/2 (due S) when η=0 at local solar noon
    ζ = mod(
        FT(3π / 2) - atan(sin(η), cos(η) * sin(ϕ) - tan(δ) * cos(ϕ)),
        FT(2π),
    )

    return d, θ, ζ
end

"""
    daily_distance_zenith_angle(
        date::DateTime,
        latitude::FT,
        (ϖ, γ, e)::Tuple{FT, FT, FT},
        param_set::IP.AIP,
    ) where {FT <: Real}

Calculates the effective zenith angle for diurnally averaged insolation and 
planet-star distance.

Returns the effective zenith angle corresponding to the diurnally
averaged insolation (averaging cos(zenith angle) over 24 hours) and 
the planet-star distance for a given date and latitude.

# Arguments
- `date::DateTime`: Current date
- `latitude::FT`: Latitude [degrees]
- `(ϖ, γ, e)::Tuple{FT, FT, FT}`: Orbital parameters tuple containing:
  - `ϖ`: Longitude of perihelion [radians]
  - `γ`: Obliquity (axial tilt) [radians]
  - `e`: Orbital eccentricity [unitless]
- `param_set::IP.AIP`: Parameter struct

# Returns
- `daily_θ`: Effective solar zenith angle [radians]
- `d`: Planet-star distance [m]
"""
function daily_distance_zenith_angle(
    date::DateTime,
    latitude::FT,
    (ϖ, γ, e)::Tuple{FT, FT, FT},
    param_set::IP.AIP;
) where {FT}
    ϕ = deg2rad(latitude)

    Δt_years = years_since_epoch(param_set, date)

    # Get distance and declination
    d, δ, _ = distance_declination_mean_anomaly(Δt_years, (ϖ, γ, e), param_set)

    # Sunrise/sunset hour angle 
    T = tan(ϕ) * tan(δ)
    # Clamp T to valid acos range and compute
    T_clamped = clamp(T, FT(-1), FT(1))
    ηd_normal = acos(-T_clamped)
    # Use ifelse to select: polar day (π), polar night (0), or normal
    ηd = ifelse(T >= FT(1), FT(π), ifelse(T <= FT(-1), FT(0), ηd_normal))

    # Effective zenith angle to get diurnally averaged insolation
    # (i.e., averaging cosine of zenith angle)
    daily_θ = mod(
        acos(FT(1 / π) * (ηd * sin(ϕ) * sin(δ) + cos(ϕ) * cos(δ) * sin(ηd))),
        FT(2π),
    )

    return daily_θ, d
end

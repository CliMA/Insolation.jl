export insolation, solar_flux_and_cos_sza

"""
    _calculate_solar_flux(d::FT, param_set::IP.AIP) where {FT <: Real}

Calculates the solar radiative energy flux at the top of the atmosphere
(TOA) based on the Earth-Sun distance and the inverse square law.

# Arguments
- `d::FT`: Earth-Sun distance [m]
- `param_set::IP.AIP`: Struct containing `tot_solar_irrad` [W m⁻²] and `orbit_semimaj` [m]
"""
function _calculate_solar_flux(d::FT, param_set::IP.AIP) where {FT <: Real}
    S0::FT = IP.tot_solar_irrad(param_set)
    d0::FT = IP.orbit_semimaj(param_set)

    # Solar radiative energy flux 
    S = S0 * (d0 / d)^2
    return S
end

"""
    insolation(θ::FT, d::FT, param_set::IP.AIP) where {FT <: Real}

Returns the top-of-atmosphere (TOA) insolation given the solar zenith
angle `θ` and Earth-Sun distance `d`.

Implements $F = S \cos(\theta)$. Insolation is set to 0 at night 
(when $\cos(\theta) < 0$).

# Arguments
- `θ::FT`: Solar zenith angle [radians]
- `d::FT`: Earth-Sun distance [m]
- `param_set::IP.AIP`: Parameter struct
"""
function insolation(θ::FT, d::FT, param_set::IP.AIP) where {FT <: Real}
    # Calculate solar radiative energy flux (W m⁻²)
    S = _calculate_solar_flux(d, param_set)

    # TOA insolation; set to 0 at night (when Sun is below horizon) 
    F = S * max(FT(0), cos(θ))
    return F
end

"""
    solar_flux_and_cos_sza(
        date::DateTime,
        od::OrbitalData,
        longitude::FT,
        latitude::FT,
        param_set::IP.AIP,
    ) where {FT <: Real}

Returns the top-of-atmosphere (TOA) solar flux (TSI weighted by
Earth-Sun distance) and the cosine of the solar zenith angle.

This format is designed for input to radiative transfer models like RRTMGP.

# Arguments
- `date::DateTime`: Current date and time
- `od::OrbitalData`: Struct with orbital parameter splines
- `longitude::FT`: Longitude [degrees]
- `latitude::FT`: Latitude [degrees]
- `param_set::IP.AIP`: Parameter struct
"""
function solar_flux_and_cos_sza(
    date::DateTime,
    od::OrbitalData,
    longitude::FT,
    latitude::FT,
    param_set::IP.AIP,
) where {FT <: Real}
    # Get instantaneous orbital parameters
    args = (
        Insolation.helper_instantaneous_zenith_angle(date, od, param_set)...,
        longitude,
        latitude,
    )
    
    # θ = solar zenith angle, ζ = solar azimuth angle, d = earth-sun distance
    θ, ζ, d = instantaneous_zenith_angle(args...)

    # Cosine of solar zenith angle (set to 0 at night) 
    μ = max(FT(0), cos(θ))

    # TOA solar flux (W m⁻²)
    S = _calculate_solar_flux(d, param_set)

    return S, μ
end

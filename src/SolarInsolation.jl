module SolarInsolation

using Dates
using ..SolarZenithAngle
using ..OrbitalParameters

export instantaneous_insolation, daily_insolation, insolation

"""
    instantaneous_insolation(date::DateTime,
                             timezone::FT,
                             longitude::FT,
                             latitude::FT) where {FT <: Real}

returns the daily averaged insolation
at a particular latitude on the given date
"""
function instantaneous_insolation(date::DateTime,
                          timezone::FT,
                          longitude::FT,
                          latitude::FT) where {FT <: Real}
    θ, d = instantaneous_zenith_angle(date, timezone, longitude, latitude)
    F = solar_insolation() * (1 / d)^2 * θ
    return F
end

"""
    daily_insolation(date::DateTime,
                     latitude::FT) where {FT <: Real}

returns the daily averaged insolation
at a particular latitude on the given date
"""
function daily_insolation(date::DateTime,
                          latitude::FT) where {FT <: Real}
    θ, d = daily_zenith_angle(date, latitude)
    F = solar_insolation() * (1 / d)^2 * θ
    return F
end

"""
    daily_insolation(days_since_equinox::I,
                       obliquity::FT,
                       perihelion::FT,
                       eccentricity::FT,
                       latitude::FT) where {FT <: Real, I <: Int}

returns the daily averaged insolation
at a particular latitude given the days since vernal equinox (defined as March 21),
orbital obliquity, longitude of perihelion, and eccentricity
"""
function daily_insolation(days_since_equinox::I,
                            obliquity::FT,
                            perihelion::FT,
                            eccentricity::FT,
                            latitude::FT) where {FT <: Real, I <: Int}
    θ, d = daily_zenith_angle(days_since_equinox, obliquity, perihelion, eccentricity, latitude)
    F = solar_insolation() * (1 / d)^2 * θ
    return F
end

"""
    insolation(θ::FT,
                     d::FT) where {FT <: Real}

returns the insolation given the zenith angle and earth-sun distance
"""
function insolation(θ::FT,
                    d::FT) where {FT <: Real}
    F = solar_insolation() * (1 / d)^2 * θ
    return F
end

end
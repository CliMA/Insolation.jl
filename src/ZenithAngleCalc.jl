module ZenithAngleCalc

using Dates
using ..OrbitalParameters

export instantaneous_zenith_angle, daily_zenith_angle

function distance_declination(date::DateTime, γ::FT, ϖ::FT, e::FT) where {FT <: Real}
    # mean anomaly, radians (3.6)
    julian_day_abs = datetime2julian(date)
    sec_since_epoch = (julian_day_abs - epoch()) * day_length()
    MA = mod(2π * sec_since_epoch / year_anom() + mean_anom_epoch(), 2π)

    # true anomaly, radians (3.8)
    TA = mod(MA + (2*e - e^FT(3/4))*sin(MA) + FT(5/4)*e^2*sin(2*MA) + FT(13/12)*e^3*sin(3*MA), 2π)

    # true longitude, radians (3.9)
    TL = mod(TA + ϖ, 2π)

    # declination, radians (3.16)
    δ = mod(asin(sin(γ) * sin(TL)), 2π)

    # earth-sun distance, (3.1)
    d = astro_unit() * (1 - e^2) / (1 + e*cos(TA))

    return d, δ
end

"""
    instantaneous_zenith_angle(date::DateTime,
                               longitude::FT,
                               latitude::FT,
                               γ::FT=obliquity(),
                               ϖ::FT=perihelion(),
                               e::FT=eccentricity()) where {FT <: Real}

returns the zenith angle and earth-sun distance
at a particular longitude and latitude on the given date (and time UTC)
given orbital parameters: obliquity, longitude of perihelion, and eccentricity
"""
function instantaneous_zenith_angle(date::DateTime,
                                    longitude::FT,
                                    latitude::FT,
                                    γ::FT=obliquity(),
                                    ϖ::FT=perihelion(),
                                    e::FT=eccentricity()) where {FT <: Real}
    λ = deg2rad(longitude)
    ϕ = deg2rad(latitude)

    d, δ = distance_declination(date, γ, ϖ, e)

    # hour angle, radians (3.17)
    julian_day_abs = datetime2julian(date)
    η_UTC = 2π * (mod(julian_day_abs, 1) - FT(0.5))
    η = η_UTC + λ

    # zenith angle, radians (3.18)
    sza = acos(cos(ϕ)*cos(δ)*cos(η) + sin(ϕ)*sin(δ))
    if sza > π/2.0
        sza = π/2.0
    end

    # solar azimuth angle
    azi = atan(sin(η), cos(η)*sin(ϕ) - tan(δ)*cos(ϕ))

    return sza, azi, d
end

"""
    daily_zenith_angle(date::DateTime,
                       latitude::FT,
                       γ::FT=obliquity(),
                       ϖ::FT=perihelion(),
                       e::FT=eccentricity()) where {FT <: Real}
returns the daily averaged zenith angle and earth-sun distance
at a particular latitude given the date and orbital parameters
obliquity, longitude of perihelion, and eccentricity
"""
function daily_zenith_angle(date::DateTime,
                            latitude::FT,
                            γ::FT=obliquity(),
                            ϖ::FT=perihelion(),
                            e::FT=eccentricity()) where {FT <: Real}
    ϕ = deg2rad(latitude)

    d, δ = distance_declination(date, γ, ϖ, e)
    
    # sunrise/sunset angle (3.19)
    T = tan(ϕ) * tan(δ)
    if T >= 1
        ηd = π
    elseif T <= -1
        ηd = 0.0
    else
        ηd = acos(-1*T)
    end
    
    # daily averaged zenith angle (3.20)
    daily_sza = acos((1/π)*(ηd*sin(ϕ)*sin(δ) + cos(ϕ)*cos(δ)*sin(ηd)))

    return daily_sza, d
end

end
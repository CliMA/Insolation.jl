module ZenithAngleCalc

using Dates
using ..OrbitalParameters

export instantaneous_zenith_angle, daily_zenith_angle

"""
    instantaneous_zenith_angle(date::DateTime,
                               timezone::FT,
                               longitude::FT,
                               latitude::FT,
                               obliquity::FT,
                               perihelion::FT,
                               eccentricity::FT) where {FT <: Real}

returns the zenith angle and earth-sun distance
at a particular longitude and latitude on the given date (and time UTC)
given orbital parameters: obliquity, longitude of perihelion, and eccentricity
"""
function instantaneous_zenith_angle(date::DateTime,
                                    longitude::FT,
                                    latitude::FT,
                                    obliquity::FT,
                                    perihelion::FT,
                                    eccentricity::FT) where {FT <: Real}
    λ = deg2rad(longitude)
    ϕ = deg2rad(latitude)
    γ = deg2rad(obliquity)
    ϖ = deg2rad(perihelion)
    e = eccentricity

    # mean anomaly at vernal equinox, radians (3.11)
    β = (1-e^2)^(0.5)
    MAve = -ϖ + (e + e^3/4)*(1+β)*sin(ϖ) - 0.5*e^2*(0.5+β)*sin(2*ϖ) + 0.25*e^3*(1/3+β)*sin(3*ϖ)

    # mean anomaly, radians (3.10)
    julian_day_abs = datetime2julian(date)
    julian_day_ref = 2451545.0 # Jan 1, 2000 at 12hr
    julian_day_ve = julian_day_ref + 79.0 # vernal equinox 2000
    t_since_vernal_equinox = mod(julian_day_abs - julian_day_ve, year_anom())
    MA = 2π * t_since_vernal_equinox / year_anom() + MAve

    # true anomaly, radians (3.8)
    TA = MA + (2*e - e^3/4)*sin(MA) + (5/4)*e^2*sin(2*MA) + (13/12)*e^3*sin(3*MA)

    # true longitude, radians (3.9)
    TL = mod(TA + ϖ, 2*π)

    # declination, radians (3.16)
    δ = mod(asin(sin(γ) * sin(TL)), 2*π)

    # hour angle, radians
    julian_day_abs = datetime2julian(date)
    η_UTC = 2π * (mod(julian_day_abs, 1.0) - 0.5)
    η = η_UTC + λ

    # zenith angle, radians (3.18)
    sza = acos(cos(ϕ)*cos(δ)*cos(η) + sin(ϕ)*sin(δ))

    # Set to a max of 90 deg.
    if sza > π/2.0
        sza = π/2.0
    end

    # solar azimuth angle
    azi = atan(sin(η), cos(η)*sin(ϕ) - tan(δ)*cos(ϕ))

    # earth-sun distance, m (3.1)
    d = astro_unit() * (1.0 - e^2) / (1.0 + e*cos(TA))

    return sza, azi, d
end

"""
    daily_zenith_angle(date::DateTime,
                       latitude::FT,
                       obliquity::FT,
                       perihelion::FT,
                       eccentricity::FT) where {FT <: Real}
returns the daily averaged zenith angle and earth-sun distance
at a particular latitude given the date, orbital obliquity, 
longitude of perihelion, and eccentricity
"""
function daily_zenith_angle(date::DateTime,
                            latitude::FT,
                            obliquity::FT,
                            perihelion::FT,
                            eccentricity::FT) where {FT <: Real}
    ϕ = deg2rad(latitude)
    γ = deg2rad(obliquity)
    ϖ = deg2rad(perihelion)
    ecc = eccentricity

    # mean anomaly at vernal equinox, radians (3.11)
    β = (1-e^2)^(0.5)
    MAve = -ϖ + (e + e^3/4)*(1+β)*sin(ϖ) - 0.5*e^2*(0.5+β)*sin(2*ϖ) + 0.25*e^3*(1/3+β)*sin(3*ϖ)

    # mean anomaly, radians (3.10)
    julian_day_abs = datetime2julian(date)
    julian_day_ref = 2451545.0 # Jan 1, 2000 at 12hr
    julian_day_ve = julian_day_ref + 79.0 # vernal equinox 2000
    t_since_vernal_equinox = mod(julian_day_abs - julian_day_ve, year_anom())
    MA = 2π * t_since_vernal_equinox / year_anom() + MAve

    # true anomaly, radians (3.8)
    TA = MA + (2*e - e^3/4)*sin(MA) + (5/4)*e^2*sin(2*MA) + (13/12)*e^3*sin(3*MA)

    # true longitude, radians (3.9)
    TL = mod(TA + ϖ, 2*π)

    # declination, radians (3.16)
    δ = mod(asin(sin(γ) * sin(TL)), 2*π)
    
    # sunrise/sunset angle
    T = tan(ϕ) * tan(δ)
    if T >= 1
        ηd = π
    elseif T <= -1
        ηd = 0.0
    else
        ηd = acos(-1*T)
    end
    
    # daily averaged zenith angle
    szabar = acos((1/π)*(ηd*sin(ϕ)*sin(δ) + cos(ϕ)*cos(δ)*sin(ηd)))

    # earth-sun distance, m (3.1)
    d = astro_unit() * (1.0 - e^2) / (1.0 + e*cos(TA))

    return szabar, d
end

end
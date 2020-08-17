module ZenithAngleCalc

using Dates
using ..OrbitalParameters

export instantaneous_zenith_angle, daily_zenith_angle

"""
    julian_century(date::DateTime)
returns the julian century (centuries since Jan 1, 2000)
given the datetime
"""
function julian_century(date::DateTime)
    # julian day
    julian_day_abs = datetime2julian(date)
    julian_day_ref = 2451545.0
    # elapsed days referenced to noon 1 Jan 2000 UTC
    jd = julian_day_abs - julian_day_ref
    # julian century
    jc = jd / (year_anom() / day_length() * 100.0)
    return jc
end

"""
    true_longitude(date::DateTime)
returns the true solar longitude in radians
given the datetime

formula from "Astronomical Algorithms" by Jean Meeus
chapter 25, ML = equation 25.2
"""
function true_longitude(date::DateTime)
    jc = julian_century(date)
    # mean solar longitude
    ML = deg2rad(mod(280.46646 + 36000.76983*jc + 0.0003032*jc^2, 360.0))
    # mean anomaly, radians
    MA = deg2rad(mod(357.52911 + 35999.05029*jc - 0.0001537*jc^2, 360.0))
    # solar equation of center
    SC = deg2rad(sin(MA)*(1.914602-0.004817*jc-0.000014*jc^2) + sin(2*MA)*(0.019993-0.000101*jc) + sin(3*MA)*0.000289)
    # true longitude
    TL = ML + SC
    return TL
end

"""
    true_anomaly(date::DateTime)
returns the true solar anomaly in radians

formula from "Astronomical Algorithms" by Jean Meeus
chapter 25, MA = equation 25.3
"""
function true_anomaly(date::DateTime)
    jc = julian_century(date)
    # mean anomaly, radians
    MA = deg2rad(mod(357.52911 + 35999.05029*jc - 0.0001537*jc^2, 360.0))
    # solar equation of center
    SC = deg2rad(sin(MA)*(1.914602-0.004817*jc-0.000014*jc^2) + sin(2*MA)*(0.019993-0.000101*jc) + sin(3*MA)*0.000289)
    # true anomaly
    TA = MA + SC
    return TA
end

"""
    eccentricity(date::DateTime)
returns the eccentricity of Earth's orbit
given the datetime

formula from "Astronomical Algorithms" by Jean Meeus
chapter 25, equation 25.4
"""
function eccentricity(date::DateTime)
    jc = julian_century(date)
    # eccentricity
    ecc = 0.016708634 - 0.000042037*jc - 0.0000001267*jc^2
    return ecc
end

"""
    obliquity(date::DateTime)
returns the obliquity of Earth's orbit in radians
given the datetime

formula from "Astronomical Algorithms" by Jean Meeus
chapter 22, approximation ignorning nutation of obliquity, Δϵ
"""
function obliquity(date::DateTime)
    jc = julian_century(date)
    # obliquity
    γ = deg2rad(mod(23.439291 - 0.01300417*jc - 1.638889e-7*jc^2 + 5.036111e-7*jc^3, 360.0))
    return γ
end

"""
    GMST(date::DateTime, timezone::FT) where {FT <: Real}

returns the Greenwich mean sidereal time in radians
given the datetime and timezone

formula from "Astronomical Algorithms" by Jean Meeus
chapter 12, equation 12.4
"""
function GMST(date::DateTime, timezone::FT) where {FT <: Real}
    jc = julian_century(date)
    # Greenwich mean sidereal time, radians
    UTC_hours = Dates.hour(date) + Dates.minute(date)/60.0 + Dates.second(date)/3600.0 - timezone
    GMST = mod(6.6974243242 + 2400.117188*jc + UTC_hours, 24.0)
    GMSTrad = mod(deg2rad(GMST * (360.0/24.0)), 2*π)
    return GMSTrad
end

"""
    earth_sun_dist(ecc::FT, TA::FT) where {FT <: Real}

returns the Earth-sun distance in meters
given the eccentricity and true anomaly

formula from "Astronomical Algorithms" by Jean Meeus
chapter 25, equation 25.5
"""
function earth_sun_dist(ecc::FT, TA::FT) where {FT <: Real}
    d_au = (1.000001018 * (1.0 - ecc^2)) / (1.0 + ecc*cos(TA))
    d = d_au * astro_unit()
    return d
end

"""
    instantaneous_zenith_angle(date::DateTime,
                               timezone::FT,
                               longitude::FT,
                               latitude::FT) where {FT <: Real}

returns the zenith angle and earth-sun distance
at a particular longitude and latitude on the given date

equations from "Astronomical Algorithms" by Jean Meeus
see documentation for details
"""
function instantaneous_zenith_angle(date::DateTime,
                                    timezone::FT,
                                    longitude::FT,
                                    latitude::FT) where {FT <: Real}
    λ = deg2rad(longitude)
    ϕ = deg2rad(latitude)
    
    TL = true_longitude(date)
    TA = true_anomaly(date)
    ecc = eccentricity(date)
    γ = obliquity(date)
    GMST = GMST(date, timezone)
    d = earth_sun_dist(ecc, TA)

    # declination, radians
    δ = mod(asin(sin(γ) * sin(TL)), 2*π)
    # right acension, radians
    RA = mod(atan(cos(γ) * sin(TL) / cos(TL)), 2*π)

    # hour angle, radians
    LMST = GMST + λ
    η = mod(LMST - RA, 2*π)

    # zenith angle
    sza = acos(cos(ϕ)*cos(δ)*cos(η) + sin(ϕ)*sin(δ))

    # parallax correction
    parallax = (planet_radius() / astro_unit()) * sin(sza)
    sza = mod(sza + parallax, 2*π)

    # Set to a max of 90 deg.
    if sza > π/2.0
        sza = π/2.0
    end

    return sza, d
end

"""
    instantaneous_zenith_angle(date::DateTime,
                               timezone::FT,
                               longitude::FT,
                               latitude::FT,
                               obliquity::FT,
                               perihelion::FT,
                               eccentricity::FT) where {FT <: Real}

returns the zenith angle and earth-sun distance
at a particular longitude and latitude on the given date
given orbital parameters: obliquity, longitude of perihelion, and eccentricity
"""
function instantaneous_zenith_angle(date::DateTime,
                                    timezone::FT,
                                    longitude::FT,
                                    latitude::FT,
                                    obliquity::FT,
                                    perihelion::FT,
                                    eccentricity::FT) where {FT <: Real}
    λ = deg2rad(longitude)
    ϕ = deg2rad(latitude)
    γ = deg2rad(obliquity)
    ϖ = deg2rad(perihelion)
    ecc = eccentricity

    TA = true_anomaly(date)
    TL = mod(TA + ϖ, 2*π)
    GMST = GMST(date, timezone)
    d = earth_sun_dist(ecc, TA)

    # declination, radians
    δ = mod(asin(sin(γ) * sin(TL)), 2*π)
    # right acension, radians
    RA = mod(atan(cos(γ) * sin(TL) / cos(TL)), 2*π)

    # hour angle, radians
    LMST = GMST + λ
    η = mod(LMST - RA, 2*π)

    # zenith angle
    sza = acos(cos(ϕ)*cos(δ)*cos(η) + sin(ϕ)*sin(δ))

    # Set to a max of 90 deg.
    if sza > π/2.0
        sza = π/2.0
    end

    return sza, d
end

"""
    daily_zenith_angle(date::DateTime,
                       latitude::FT) where {FT <: Real}

returns the daily averaged zenith angle and earth-sun distance
at a particular latitude on the given date
"""
function daily_zenith_angle(date::DateTime,
                            latitude::FT) where {FT <: Real}
    ϕ = deg2rad(latitude)

    TL = true_longitude(date)
    TA = true_anomaly(date)
    ecc = eccentricity(date)
    γ = obliquity(date)
    d = earth_sun_dist(ecc, TA)

    # declination, radians
    δ = mod(asin(sin(γ) * sin(TL)), 2*π)
    
    # sunrise/sunset angle
    T = tan(ϕ) * tan(δ)
    if T >= 1
        ηd = π
    elseif T <= -1
        ηd = 0.0
    else
        ηd = acos(-1 * T)
    end
    
    # daily averaged zenith angle
    szabar = acos((1/π)*(ηd*sin(ϕ)*sin(δ) + cos(ϕ)*cos(δ)*sin(ηd)))

    return szabar, d
end

"""
    daily_zenith_angle(days_since_equinox::I,
                       obliquity::FT,
                       perihelion::FT,
                       eccentricity::FT,
                       latitude::FT) where {FT <: Real, I <: Int}

returns the daily averaged zenith angle and earth-sun distance
at a particular latitude given the days since vernal equinox (defined as March 21),
orbital obliquity, longitude of perihelion, and eccentricity
"""
function daily_zenith_angle(days_since_equinox::I,
                            obliquity::FT,
                            perihelion::FT,
                            eccentricity::FT,
                            latitude::FT) where {FT <: Real, I <: Int}
    γ = deg2rad(obliquity)
    ϖ = deg2rad(perihelion)
    ecc = eccentricity
    ϕ = deg2rad(latitude)

    # calculate mean anomaly at vernal equinox and mean anomaly
    β = (1 - ecc^2)^0.5
    MA_VE = mod(-ϖ + (ecc + ecc^3/4)*(1+β)*sin(ϖ), 2*π)
    MA = mod(2*π*days_since_equinox / (year_anom() / day_length()) + MA_VE, 2*π)

    # solar longitude and true anomaly
    TA = mod(MA + (2*ecc - ecc^3/4)*sin(MA), 2*π)
    TL = mod(TA + ϖ, 2*π)

    # radius earth-sun distance, AU and m
    d = earth_sun_dist(ecc, TA)

    # declination, radians
    δ = mod(asin(sin(γ) * sin(TL)), 2*π)
    
    # sunrise/sunset angle
    T = tan(ϕ) * tan(δ)
    if T >= 1
        ηd = π
    elseif T <= -1
        ηd = 0.0
    else
        ηd = acos(-1 * T)
    end
    
    # daily averaged zenith angle
    szabar = acos((1/π)*(ηd*sin(ϕ)*sin(δ) + cos(ϕ)*cos(δ)*sin(ηd)))

    return szabar, d
end

end
module SolarZenithAngle

using Dates
using ..OrbitalParameters

export instantaneous_zenith_angle, daily_zenith_angle

"""
    instantaneous_zenith_angle(date::DateTime,
                               timezone::FT,
                               longitude::FT,
                               latitude::FT) where {FT <: Real}

returns the zenith angle and earth-sun distance
at a particular longitude and latitude on the given date

add citations: 
- https://www.esrl.noaa.gov/gmd/grad/solcalc/
- https://www.cfa.harvard.edu/~jzhao/times.html
- https://github.com/thabbott/zenithangle/blob/master/solar.js
- https://github.com/claresinger/3d-cloud-rad/blob/master/sza-scripts/sza_utils.py
- http://farside.ph.utexas.edu/Books/Syntaxis/Almagest/node36.html
- https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
- "Astronomical Algorithms" by Jean Meeus
"""
function instantaneous_zenith_angle(date::DateTime,
                                    timezone::FT,
                                    longitude::FT,
                                    latitude::FT) where {FT <: Real}
    λ = deg2rad(longitude)
    ϕ = deg2rad(latitude)
    
    # julian day
    julian_day_abs = datetime2julian(date)
    # elapsed days referenced to noon 1 Jan 2000 UTC
    jd = julian_day_abs - 2451545.0
    # julian century
    jc = jd / (year_anom() / day_length() * 100.0)

    # mean solar longitude, radians
    ML = deg2rad(mod(280.46646 + 36000.76983*jc + 0.0003032*jc^2, 360.0))
    # mean anomaly, radians
    MA = deg2rad(mod(357.52911 + 35999.05029*jc - 0.0001537*jc^2, 360.0))
    # solar equation of center
    SC = deg2rad(sin(MA)*(1.914602-0.004817*jc-0.000014*jc^2) + sin(2*MA)*(0.019993-0.000101*jc) + sin(3*MA)*0.000289)
    # true longitude
    TL = ML + SC
    # true anomaly
    TA = MA + SC

    # eccentricity
    ecc = 0.016708634 - 0.000042037*jc - 0.0000001267*jc^2
    # obliquity, radians
    γ = deg2rad(mod(23.439291 - 0.01300417*jc - 1.638889e-7*jc^2 + 5.036111e-7*jc^3, 360.0))
    # longitude of perihelion, radians
    ϖ = mod(TL - TA, 2*π)

    # radius earth-sun distance, AU and m
    d_au = (1.000001018 * (1.0 - ecc^2)) / (1.0 + ecc*cos(TA))
    d = d_au * astro_unit()

    # declination, radians
    δ = mod(asin(sin(γ) * sin(TL)), 2*π)
    # right acension, radians
    RA = mod(atan(cos(γ) * sin(TL) / cos(TL)), 2*π)

    # hour angle, radians
    UTC_hours = Dates.hour(date) + Dates.minute(date)/60.0 + Dates.second(date)/3600.0 - timezone
    GMST = mod(6.6974243242 + 2400.117188*jc + UTC_hours, 24.0)
    GMSTrad = mod(deg2rad(GMST * (360.0/24.0)), 2*π)
    LMST = GMSTrad + λ
    η = mod(LMST - RA, 2*π)

    # zenith angle
    sza = acos(cos(ϕ)*cos(δ)*cos(η) + sin(ϕ)*sin(δ))

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

    # julian day
    julian_day_abs = datetime2julian(date)
    # elapsed days referenced to noon 1 Jan 2000 UTC
    jd = julian_day_abs - 2451545.0
    # julian century
    jc = jd / (year_anom() / day_length() * 100.0)
    
    # days since equinox
    days_since_equinox = Dates.day(date) - 76.0

    # calculate mean anomaly at vernal equinox and mean anomaly
    β = (1 - ecc^2)^0.5
    MA_VE = mod(-ϖ + (ecc + ecc^3/4)*(1+β)*sin(ϖ), 2*π)
    MA = mod(2*π*days_since_equinox / (year_anom() / day_length()) + MA_VE, 2*π)

    # solar longitude and true anomaly
    TA = mod(MA + (2*ecc - ecc^3/4)*sin(MA), 2*π)
    TL = mod(TA + ϖ, 2*π)

    # radius earth-sun distance, AU and m
    d_au = (1.000001018 * (1.0 - ecc^2)) / (1.0 + ecc*cos(TA))
    d = d_au * astro_unit()

    # declination, radians
    δ = mod(asin(sin(γ) * sin(TL)), 2*π)
    # right acension, radians
    RA = mod(atan(cos(γ) * sin(TL) / cos(TL)), 2*π)

    # hour angle, radians
    UTC_hours = Dates.hour(date) + Dates.minute(date)/60.0 + Dates.second(date)/3600.0 - timezone
    GMST = mod(6.6974243242 + 2400.117188*jc + UTC_hours, 24.0)
    GMSTrad = mod(deg2rad(GMST * (360.0/24.0)), 2*π)
    LMST = GMSTrad + λ
    η = mod(LMST - RA, 2*π)

    # zenith angle
    sza = acos(cos(ϕ)*cos(δ)*cos(η) + sin(ϕ)*sin(δ))

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
    
    # julian day
    julian_day_abs = datetime2julian(date)
    # elapsed days referenced to noon 1 Jan 2000 UTC
    jd = julian_day_abs - 2451545.0
    # julian century
    jc = jd / (year_anom() / day_length() * 100.0)

    # mean solar longitude, radians
    ML = deg2rad(mod(280.46646 + 36000.76983*jc + 0.0003032*jc^2, 360.0))
    # mean anomaly, radians
    MA = deg2rad(mod(357.52911 + 35999.05029*jc - 0.0001537*jc^2, 360.0))
    # solar equation of center
    SC = deg2rad(sin(MA)*(1.914602-0.004817*jc-0.000014*jc^2) + sin(2*MA)*(0.019993-0.000101*jc) + sin(3*MA)*0.000289)
    # true longitude
    TL = ML + SC
    # true anomaly
    TA = MA + SC

    # eccentricity
    ecc = 0.016708634 - 0.000042037*jc - 0.0000001267*jc^2
    # obliquity, radians
    γ = deg2rad(mod(23.439291 - 0.01300417*jc - 1.638889e-7*jc^2 + 5.036111e-7*jc^3, 360.0))
    # longitude of perihelion, radians
    ϖ = mod(TL - TA, 2*π)

    # radius earth-sun distance, AU and m
    d_au = (1.000001018 * (1.0 - ecc^2)) / (1.0 + ecc*cos(TA))
    d = d_au * astro_unit()

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
    d_au = (1.000001018 * (1.0 - ecc^2)) / (1.0 + ecc*cos(TA))
    d = d_au * astro_unit()

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
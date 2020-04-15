module SZA

using Dates
include("OrbitalParameters.jl")

export instantaneous_zenith_angle, daily_zenith_angle

"""
    instantaneous_zenith_angle(date::DateTime,
                               timezone::FT,
                               longitude::FT,
                               latitude::FT) where {FT <: Real}

    returns the zenith angle at a particular longitude and latitude on the given date
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
    jc = jd / (OrbitalParameters.year_anom / OrbitalParameters.day_length * 100.0)

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
    d = d_au * OrbitalParameters.astro_unit

    # declination, radians
    δ = mod(asin(sin(γ) * sin(TL)), 2*π)
    # right acension, radians
    RA = mod(atan(cos(γ) * sin(TL) / cos(TL)), 2*π)

    # hour angle, radians
    hours = Dates.hour(date) - timezone
    GMST = mod(6.6974243242 + 2400.117188*jc + hours, 24.0)
    GMSTrad = mod(deg2rad(GMST * (360.0/24.0)), 2*π)
    LMST = GMSTrad + λ
    η = mod(LMST - RA, 2*π)

    # zenith angle
    sza = acos(cos(ϕ)*cos(δ)*cos(η) + sin(ϕ)*sin(δ))

    return sza
end

"""
    daily_zenith_angle(date::DateTime,
                       longitude::FT,
                       latitude::FT) where {FT <: Real}

    returns the daily averaged zenith angle at a particular longitude and latitude on the given date
"""
function daily_zenith_angle(date::DateTime,
                      longitude::FT,
                      latitude::FT) where {FT <: Real}
    λ = deg2rad(longitude)
    ϕ = deg2rad(latitude)
    
    # julian day
    julian_day_abs = datetime2julian(date)
    # elapsed days referenced to noon 1 Jan 2000 UTC
    jd = julian_day_abs - 2451545.0
    # julian century
    jc = jd / (OrbitalParameters.year_anom / OrbitalParameters.day_length * 100.0)

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
    d = d_au * OrbitalParameters.astro_unit

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

    return szabar
end


end
module SZA
using Dates
include("OrbitalParameters.jl")

"""
    zenith angle(date::DateTime,
                longitude::FT,
                latitude::FT) where {FT <: Real}

    returns the zenith angle at a particular longitude and latitude on the given date
"""
function zenith_angle(date::DateTime,
                      longitude::FT,
                      latitude::FT) where {FT <: Real}
    lonrad = deg2rad(longitude)
    latrad = deg2rad(latitude)
    
    julian_day = datetime2julian(date)
    elapsed_julian_days = julian_day - 2451545.0 # elapsed since noon 1 Jan 2000 UTC
    hours = Dates.hour(date) + (Dates.minute(date) + Dates.second(date) / 60.) / 60.

    # orbital parameters for Earth in present day
    #γ = 0.4090928 - 6.2140e-9 * elapsed_julian_days + 0.0000396 * cos(2.1429 - 0.0010394594 * elapsed_julian_days)
    γ = deg2rad(23.44)
    ϖ = deg2rad(282.95)
    ecc = 0.0167086

    # celestial coordinates, declination and hour angle
    # dMeanLongitude = mod(4.8950630 + 0.017202791698 * elapsed_julian_days, 2*π)  # Radians
    # dMeanAnomaly = mod(6.2400600 + 0.0172019699 * elapsed_julian_days, 2*π) # Radians
    # dOmega = mod(2.1429 - 0.0010394594 * elapsed_julian_days, 2*π) # Radians
    # dEclipticLongitude = dMeanLongitude + (2*ecc - ecc^3/4) * sin(dMeanAnomaly)\
    #     + (5*ecc^2/4) * sin(2. * dMeanAnomaly) - 0.0001134 - 0.0000203 * sin(dOmega)

    dOmega = mod(2.1429 - 0.0010394594 * elapsed_julian_days, 2*π) # Radians
    dϖ = mod(4.8950630 + 0.017202791698 * elapsed_julian_days - 0.0001134 - 0.0000203 * sin(dOmega), 2*π)
    dM = mod(6.2400600 + 0.0172019699 * elapsed_julian_days, 2*π)
    dA = mod(dM + (2*ecc - ecc^3/4) * sin(dM) + (5*ecc^2/4) * sin(2*dM), 2*π)
    dLs = mod(dϖ + dA, 2*π)

    M_VE = mod(- ϖ + (ecc + ecc^3 / 4) * (1 + sqrt(1 - ecc^2)) * sin(ϖ), 2*π)
    t = Dates.day(date) + (Dates.month(date)-1)*30.0
    M = mod((2 * π * (t - 76.0)) / (OrbitalParameters.year_anom) + M_VE, 2*π)
    A = mod(M + (2*ecc - ecc^3/4) * sin(M) + (5*ecc^2/4) * sin(2*M), 2*π)
    Ls = mod(ϖ + A, 2*π)

    println(dϖ)
    println(ϖ)
    println()

    println(dM)
    println(M)
    println()

    println(dA)
    println(A)
    println()

    println(dLs)
    println(Ls)
    println()

    # declination = arcsin(sin(γ) * sin(Ls))


    # right_acension = 
    # greenwich_sidereal_time = 6.6974243242 + 0.0657098283 * elapsed_julian_days + hours
    # local_sidereal_time = rad2deg(greenwich_sidereal_time * 15.0 + longitude)
    # hour_angle = local_sidereal_time - right_acension

    # sza = (arccos(cos(latrad) * cos(hour_angle) * cos(declination) + sin(declination) * sin(latrad)))
    # parallax = (OrbitalParameters.planet_radius / OrbitalParameters.astro_unit) * sin(sza)
    # sza = rad2deg(sza + parallax)

    sza = 9.0
    return sza
end

function main()
    rightnow = Dates.now()
    sza = zenith_angle(rightnow, -120.0, 30.0)
    println(sza)
end

main()

end

# γ =                 deg2rad(23.44)          # rad
# ϖ =                 deg2rad(282.95)         # rad
# e =                 0.016708035             # None
export instantaneous_zenith_angle, daily_zenith_angle

function mean_anomaly_vernal_equinox(ϖ::FT, e::FT) where {FT <: Real}
    # mean anomaly at vernal equinox, radians (3.11)
    β = (1-e^2)^0.5
    M_v = -ϖ + (e+0.25*e^3)*(1+β)*sin(ϖ) - 0.5*e^2*(0.5+β)*sin(2ϖ) + 0.25*e^3*(1/3+β)*sin(3ϖ)
    return M_v
end

function distance_declination(date::DateTime, param_set::APS, γ::FT, ϖ::FT, e::FT) where {FT <: Real}
    Ya::FT = year_anom(param_set)
    day::FT = Planet.day(param_set)
    AU::FT = astro_unit()
    
    # time of vernal equinox in the epoch (rearrangement of 3.6 and 3.10)
    M_v0 = mean_anomaly_vernal_equinox(ϖ_epoch(), e)
    time_v = Ya * (M_v0 - M_epoch()) / 2π + epoch()*day

    # mean anomaly given mean anomaly at vernal equinox (3.10)
    time = datetime2julian(date)*day
    M_v = mean_anomaly_vernal_equinox(ϖ, e)
    MA = mod(2π * (time - time_v) / Ya + M_v, 2π)

    # true anomaly, radians (3.8)
    TA = mod(MA + (2*e - e^FT(3/4))*sin(MA) + FT(5/4)*e^2*sin(2*MA) + FT(13/12)*e^3*sin(3*MA), 2π)

    # true longitude, radians (3.9)
    TL = mod(TA + ϖ, 2π)

    # declination, radians (3.16)
    δ = mod(asin(sin(γ) * sin(TL)), 2π)

    # earth-sun distance, (3.1)
    d = AU * (1 - e^2) / (1 + e*cos(TA))

    return d, δ
end

"""
    instantaneous_zenith_angle(date::DateTime,
                               longitude::FT,
                               latitude::FT,
                               param_set::APS,
                               γ::FT=γ_epoch(),
                               ϖ::FT=ϖ_epoch(),
                               e::FT=e_epoch()) where {FT <: Real}

returns the zenith angle and earth-sun distance
at a particular longitude and latitude on the given date (and time UTC)
given orbital parameters: obliquity, longitude of perihelion, and eccentricity
param_set is an AbstractParameterSet from CLIMAParameters.jl
"""
function instantaneous_zenith_angle(date::DateTime,
                                    longitude::FT,
                                    latitude::FT,
                                    param_set::APS,
                                    γ::FT=γ_epoch(),
                                    ϖ::FT=ϖ_epoch(),
                                    e::FT=e_epoch()) where {FT <: Real}
    λ = deg2rad(longitude)
    ϕ = deg2rad(latitude)

    d, δ = distance_declination(date, param_set, γ, ϖ, e)

    # hour angle, zero at local solar noon, radians (3.17)
    julian_day_abs = datetime2julian(date)
    η_UTC = 2π * (mod(julian_day_abs, 1) - FT(0.5))
    η = mod(η_UTC + λ, 2π)

    # zenith angle, radians (3.18)
    sza = mod(acos(cos(ϕ)*cos(δ)*cos(η) + sin(ϕ)*sin(δ)), 2π)

    # solar azimuth angle, azi = 0 when due E and increasing CCW
    # azi = 3π/2 (due S) when η=0 at local solar noon
    azi = mod(3π/2 - atan(sin(η), cos(η)*sin(ϕ) - tan(δ)*cos(ϕ)), 2π)

    return sza, azi, d
end

"""
    daily_zenith_angle(date::DateTime,
                       latitude::FT,
                       param_set::APS,
                       γ::FT=γ_epoch(),
                       ϖ::FT=ϖ_epoch(),
                       e::FT=e_epoch()) where {FT <: Real}
returns the daily averaged zenith angle and earth-sun distance
at a particular latitude given the date and orbital parameters
obliquity, longitude of perihelion, and eccentricity
param_set is an AbstractParameterSet from CLIMAParameters.jl
"""
function daily_zenith_angle(date::DateTime,
                            latitude::FT,
                            param_set::APS,
                            γ::FT=γ_epoch(),
                            ϖ::FT=ϖ_epoch(),
                            e::FT=e_epoch()) where {FT <: Real}
    ϕ = deg2rad(latitude)

    d, δ = distance_declination(date, param_set, γ, ϖ, e)
    
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
    daily_sza = mod(acos((1/π)*(ηd*sin(ϕ)*sin(δ) + cos(ϕ)*cos(δ)*sin(ηd))), 2π)

    return daily_sza, d
end
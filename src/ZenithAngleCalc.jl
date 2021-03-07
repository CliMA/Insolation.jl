export instantaneous_zenith_angle, daily_zenith_angle

# mean anomaly at vernal equinox, radians (eq 3.11)
function mean_anomaly_vernal_equinox(ϖ::FT, e::FT) where {FT <: Real}
    β = (FT(1)-e^FT(2))^FT(1/2)
    M_v = -ϖ + (e+FT(1/4)*e^FT(3))*(FT(1)+β)*sin(ϖ)
            - FT(1/2)*e^FT(2)*(FT(1/2)+β)*sin(2ϖ)
            + FT(1/4)*e^FT(3)*(FT(1/3)+β)*sin(3ϖ)
    M_v = mod(M_v, 2π)
    return M_v
end

# true anomaly, radians (eq 3.8)
function true_anomaly(MA::FT, e::FT) where {FT <: Real}
    TA = MA + (FT(2)*e - FT(1/4)*e^FT(3))*sin(MA) 
            + FT(5/4)*e^2*sin(2MA) 
            + FT(13/12)*e^3*sin(3MA)
    TA = mod(TA, 2π)
    return TA
end

function distance_declination(::Type{FT}, date::DateTime, param_set::APS) where {FT <: Real}
    Ya::FT = year_anom(param_set)
    day_length::FT = Planet.day(param_set)
    AU::FT = astro_unit()

    _epoch::FT = epoch(param_set)
    M0::FT = mean_anom_epoch(param_set)
    ϖ0::FT = lon_perihelion_epoch(param_set)

    γ::FT = obliq_epoch(param_set)
    ϖ::FT = lon_perihelion(param_set)
    e::FT = eccentricity_epoch(param_set)
    
    # time of vernal equinox in the epoch (rearrangement of 3.6 and 3.10)
    M_v0 = mean_anomaly_vernal_equinox(ϖ0, e)
    time_v = Ya * (M_v0 - M0) / 2π + _epoch

    # mean anomaly given mean anomaly at vernal equinox (3.10)
    time = datetime2julian(date)*day_length
    M_v = mean_anomaly_vernal_equinox(ϖ, e)
    MA = mod(2π * (time - time_v) / Ya + M_v, 2π)

    # true anomaly, radians (3.8)
    TA = true_anomaly(MA, e)

    # true longitude, radians (3.9)
    TL = mod(TA + ϖ, 2π)

    # declination, radians (3.16)
    δ = mod(asin(sin(γ) * sin(TL)), 2π)

    # earth-sun distance, (3.1)
    d = AU * (1 - e^FT(2)) / (FT(1) + e*cos(TA))

    return d, δ
end

"""
    instantaneous_zenith_angle(date::DateTime,
                               longitude::FT,
                               latitude::FT,
                               param_set::APS) where {FT <: Real}

returns the zenith angle and earth-sun distance
at a particular longitude and latitude on the given date (and time UTC)
given orbital parameters: obliquity, longitude of perihelion, and eccentricity
param_set is an AbstractParameterSet from CLIMAParameters.jl
"""
function instantaneous_zenith_angle(date::DateTime,
                                    longitude::FT,
                                    latitude::FT,
                                    param_set::APS) where {FT <: Real}
    λ = deg2rad(longitude)
    ϕ = deg2rad(latitude)

    d, δ = distance_declination(FT, date, param_set)

    # hour angle, zero at local solar noon, radians (3.17)
    julian_day_abs = datetime2julian(date)
    η_UTC = 2π * mod(julian_day_abs, 1)
    η = mod(η_UTC + λ, 2π)

    # zenith angle, radians (3.18)
    θ = mod(acos(cos(ϕ)*cos(δ)*cos(η) + sin(ϕ)*sin(δ)), 2π)

    # solar azimuth angle, ζ = 0 when due E and increasing CCW
    # ζ = 3π/2 (due S) when η=0 at local solar noon
    ζ = mod(3π/2 - atan(sin(η), cos(η)*sin(ϕ) - tan(δ)*cos(ϕ)), 2π)

    return θ, ζ, d
end

"""
    daily_zenith_angle(date::DateTime,
                       latitude::FT,
                       param_set::APS) where {FT <: Real}
returns the daily averaged zenith angle and earth-sun distance
at a particular latitude given the date and orbital parameters
obliquity, longitude of perihelion, and eccentricity
param_set is an AbstractParameterSet from CLIMAParameters.jl
"""
function daily_zenith_angle(date::DateTime,
                            latitude::FT,
                            param_set::APS) where {FT <: Real}
    ϕ = deg2rad(latitude)

    d, δ = distance_declination(FT, date, param_set)
    
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
    daily_θ = mod(acos((1/π)*(ηd*sin(ϕ)*sin(δ) + cos(ϕ)*cos(δ)*sin(ηd))), 2π)

    return daily_θ, d
end
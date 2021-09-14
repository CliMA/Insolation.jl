export instantaneous_zenith_angle, daily_zenith_angle

# mean anomaly at vernal equinox, radians (eq 3.11)
function mean_anomaly_vernal_equinox(ϖ::FT, e::FT) where {FT <: Real}
    β = (FT(1)-e^FT(2))^FT(0.5)
    M_v = -ϖ + (e+FT(0.25)*e^FT(3))*(FT(1)+β)*sin(ϖ)
            - FT(0.5)*e^FT(2)*(FT(0.5)+β)*sin(FT(2)*ϖ)
            + FT(0.25)*e^FT(3)*(FT(1)/FT(3)+β)*sin(FT(3)*ϖ)
    M_v = mod(M_v, FT(2)*FT(π))
    return M_v
end

# true anomaly, radians (eq 3.8)
function true_anomaly(MA::FT, e::FT) where {FT <: Real}
    TA = MA + (FT(2)*e - FT(0.25)*e^FT(3))*sin(MA) 
            + FT(1.25)*e^FT(2)*sin(FT(2)*MA) 
            + FT(13)/FT(12)*e^FT(3)*sin(FT(3)*MA)
    TA = mod(TA, FT(2)*FT(π))
    return TA
end

# equation of time, radians
function equation_of_time(e::FT, MA::FT, γ::FT, ϖ::FT) where {FT <: Real}
    Δt = -FT(2)*e*sin(MA) + tan(γ/FT(2))^FT(2)*sin(FT(2)*(MA + ϖ)) 
    Δt = mod(Δt+π, FT(2)*FT(π))-π
end

# calculate the distance, declination, and hour angle (at lon=0)
function distance_declination_hourangle(::Type{FT}, date::DateTime, param_set::APS, eot_correction::Bool) where {FT <: Real}
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
    time_v = Ya * (M_v0 - M0) / (FT(2)*FT(π)) + _epoch

    # mean anomaly given mean anomaly at vernal equinox (3.10)
    time = datetime2julian(date) * day_length # in seconds
    M_v = mean_anomaly_vernal_equinox(ϖ, e)
    MA = mod(FT(2)*FT(π) * FT(time - time_v) / Ya + M_v, FT(2)*FT(π))

    # true anomaly, radians (3.8)
    TA = true_anomaly(MA, e)

    # true longitude, radians (3.9)
    TL = mod(TA + ϖ, FT(2)*FT(π))

    # declination, radians (3.16)
    δ = mod(asin(sin(γ) * sin(TL)), FT(2)*FT(π))

    # earth-sun distance, (3.1)
    d = AU * (FT(1) - e^FT(2)) / (FT(1) + e*cos(TA))

    # hour angle, zero at local solar noon, radians (3.17)
    if eot_correction
        Δt = equation_of_time(e, MA, γ, ϖ) / 2π * day_length # radians to seconds
    else
        Δt = FT(0)
    end
    η_UTC = mod(FT(2)*FT(π) * FT(time + Δt) / day_length, FT(2)*FT(π))

    return d, δ, η_UTC
end

"""
    instantaneous_zenith_angle(date::DateTime,
                               longitude::FT,
                               latitude::FT,
                               param_set::APS;
                               eot_correction::Bool=true) where {FT <: Real}

Returns the zenith angle and earth-sun distance
at a particular longitude and latitude on the given date (and time UTC)
given orbital parameters: obliquity, longitude of perihelion, and eccentricity
param_set is an AbstractParameterSet from CLIMAParameters.jl.

eot_correction is an optional Boolean keyword argument that defaults to true
when set to true the equation of time correction is turned on.
This switch functionality is implemented for easy comparisons with reanalyses.
"""
function instantaneous_zenith_angle(date::DateTime,
                                    longitude::FT,
                                    latitude::FT,
                                    param_set::APS; 
                                    eot_correction::Bool=true) where {FT <: Real}
    λ = deg2rad(longitude)
    ϕ = deg2rad(latitude)

    d, δ, η_UTC = distance_declination_hourangle(FT, date, param_set, eot_correction)

    # hour angle
    η = mod(η_UTC + λ, FT(2)*FT(π))

    # zenith angle, radians (3.18)
    θ = mod(acos(cos(ϕ)*cos(δ)*cos(η) + sin(ϕ)*sin(δ)), FT(2)*FT(π))

    # solar azimuth angle, ζ = 0 when due E and increasing CCW
    # ζ = 3π/2 (due S) when η=0 at local solar noon
    ζ = mod(FT(1.5)*FT(π) - atan(sin(η), cos(η)*sin(ϕ) - tan(δ)*cos(ϕ)), FT(2)*FT(π))

    return θ, ζ, d
end

"""
    daily_zenith_angle(date::DateTime,
                       latitude::FT,
                       param_set::APS;
                       eot_correction::Bool=true) where {FT <: Real}
Returns the daily averaged zenith angle and earth-sun distance
at a particular latitude given the date and orbital parameters
obliquity, longitude of perihelion, and eccentricity
param_set is an AbstractParameterSet from CLIMAParameters.jl.

eot_correction is an optional Boolean keyword argument that defaults to true
when set to true the equation of time correction is turned on.
This switch functionality is implemented for easy comparisons with reanalyses.
"""
function daily_zenith_angle(date::DateTime,
                            latitude::FT,
                            param_set::APS;
                            eot_correction::Bool=true) where {FT <: Real}
    ϕ = deg2rad(latitude)

    d, δ, _ = distance_declination_hourangle(FT, date, param_set, eot_correction)
    
    # sunrise/sunset angle (3.19)
    T = tan(ϕ) * tan(δ)
    if T >= FT(1)
        ηd = FT(π)
    elseif T <= FT(-1)
        ηd = FT(0.0)
    else
        ηd = acos(FT(-1)*T)
    end
    
    # daily averaged zenith angle (3.20)
    daily_θ = mod(acos(FT(1)/FT(π)*(ηd*sin(ϕ)*sin(δ) + cos(ϕ)*cos(δ)*sin(ηd))), FT(2)*FT(π))

    return daily_θ, d
end
export instantaneous_zenith_angle, daily_zenith_angle

# true anomaly, radians
function true_anomaly(MA::FT, e::FT) where {FT <: Real}
    TA = MA + (2 * e - FT(1 / 4) * e^3) * sin(MA)
    +FT(5 / 4) * e^2 * sin(2 * MA)
    +FT(13 / 12) * e^3 * sin(3 * MA)
    return mod(TA, FT(2π))
end

# equation of time, radians
function equation_of_time(e::FT, MA::FT, γ::FT, ϖ::FT) where {FT <: Real}
    _Δt = -2 * e * sin(MA) + tan(γ / 2)^2 * sin(2 * (MA + ϖ))
    return mod(_Δt + FT(π), FT(2π)) - FT(π)
end

# calculate orbital parameters (Milankovitch)
compute_orbital_parameters(od, Δt_years::FT) where {FT} = (
    FT(ϖ_spline(od, Δt_years)),
    FT(γ_spline(od, Δt_years)),
    FT(e_spline(od, Δt_years)))

compute_orbital_parameters(param_set) = (
    IP.lon_perihelion_epoch(param_set),
    IP.obliq_epoch(param_set),
    IP.eccentricity_epoch(param_set))

function get_Δt_years(
    param_set::IP.InsolationParameters{FT},
    date,
    epoch_string="2000-01-01T11:58:56.816",
) where {FT}
    (; year_anom, day) = param_set
    date0 = DateTime(epoch_string,dateformat"y-m-dTHH:MM:SS.s")
    days_per_year = year_anom / day
    return FT(datetime2julian(date) - datetime2julian(date0)) / days_per_year
end

# calculate the distance, declination, and hour angle (at lon=0)
function distance_declination_hourangle(date::DateTime,
                                        (ϖ, γ, e)::Tuple{FT,FT,FT},
                                        epoch_string,
                                        param_set::IP.AIP,
                                        eot_correction::Bool,
) where {FT}
    Ya = IP.year_anom(param_set)
    day_length = IP.day(param_set)
    d0 = IP.orbit_semimaj(param_set)
    M0 = IP.mean_anom_epoch(param_set)

    date0 = DateTime(epoch_string,dateformat"y-m-dTHH:MM:SS.s")

    days_per_year = Ya / day_length
    Δt_years =
        FT(datetime2julian(date) - datetime2julian(date0)) / days_per_year

    # mean anomaly given mean anomaly at epoch
    MA = mod(FT(2π) * (Δt_years) + M0, FT(2π))

    # true anomaly, radians
    TA = true_anomaly(MA, e)

    # true longitude, radians
    TL = mod(TA + ϖ, FT(2π))

    # declination, radians
    δ = mod(asin(sin(γ) * sin(TL)), FT(2π))

    # earth-sun distance
    d = d0 * (1 - e^2) / (1 + e * cos(TA))

    # hour angle, zero at local solar noon, radians
    if eot_correction
        Δt = equation_of_time(e, MA, γ, ϖ) / FT(2π) * day_length # radians to seconds
    else
        Δt = FT(0)
    end
    time_in_day = FT(mod(datetime2julian(date), 1)) * day_length # s in this day
    η_UTC = mod(FT(2π) * (time_in_day + Δt) / day_length, FT(2π))

    return d, δ, η_UTC
end

"""
    instantaneous_zenith_angle(date::DateTime,
                               od::OrbitalData,
                               longitude::FT,
                               latitude::FT,
                               param_set::IP.AIP;
                               eot_correction::Bool=true,
                               milankovitch::Bool=true) where {FT <: Real}

Returns the zenith angle and earth-sun distance
at a particular longitude and latitude on the given date (and time UTC)
given orbital parameters: obliquity, longitude of perihelion, and eccentricity
param_set is an AbstractParameterSet from CLIMAParameters.jl.

`eot_correction` is an optional Boolean keyword argument that defaults to true
when set to true the equation of time correction is turned on.
This switch functionality is implemented for easy comparisons with reanalyses.

`milankovitch` is an optional Boolean keyword argument that defaults to true
when set to true the orbital parameters are calculated for the given DateTime
when set to false the orbital parameters at the J2000 epoch from CLIMAParameters are used.
"""
function instantaneous_zenith_angle(date::DateTime,
                                    epoch_string,
                                    (ϖ, γ, e)::Tuple{FT,FT,FT},
                                    longitude::FT,
                                    latitude::FT,
                                    param_set::IP.AIP,
                                    eot_correction::Bool,
) where {FT}
    λ = deg2rad(longitude)
    ϕ = deg2rad(latitude)


    d, δ, η_UTC = distance_declination_hourangle(date, (ϖ, γ, e), epoch_string, param_set, eot_correction)

    # hour angle
    η = mod(η_UTC + λ, FT(2π))

    # zenith angle, radians
    θ = mod(acos(cos(ϕ) * cos(δ) * cos(η) + sin(ϕ) * sin(δ)), FT(2π))

    # solar azimuth angle, ζ = 0 when due E and increasing CCW
    # ζ = 3π/2 (due S) when η=0 at local solar noon
    ζ = mod(
        FT(3π / 2) - atan(sin(η), cos(η) * sin(ϕ) - tan(δ) * cos(ϕ)),
        FT(2π),
    )

    return θ, ζ, d
end

function instantaneous_zenith_angle(date::DateTime,
                                    od::OrbitalData,
                                    longitude::FT,
                                    latitude::FT,
                                    param_set::IP.AIP;
                                    eot_correction::Bool=true,
                                    ) where {FT}

    epoch_string = "2000-01-01T11:58:56.816"
    # epoch_string::String = IP.epoch(param_set)
    Δt_years = get_Δt_years(param_set, date, epoch_string)
    instantaneous_zenith_angle(
        date,
        epoch_string,
        compute_orbital_parameters(od, Δt_years),
        longitude,
        latitude,
        param_set,
        eot_correction,
    )
end

function instantaneous_zenith_angle(date::DateTime,
                                    longitude::FT,
                                    latitude::FT,
                                    param_set::IP.AIP;
                                    eot_correction::Bool=true,
                                    ) where {FT}
    epoch_string = "2000-01-01T11:58:56.816"
    # epoch_string::String = IP.epoch(param_set)

    return instantaneous_zenith_angle(
        date,
        epoch_string,
        compute_orbital_parameters(param_set),
        longitude,
        latitude,
        param_set,
        eot_correction,
    )

end

"""
    daily_zenith_angle(date::DateTime,
                       od::OrbitalData,
                       latitude::FT,
                       param_set::IP.AIP;
                       eot_correction::Bool=true,
                       milankovitch::Bool=true) where {FT <: Real}

Returns the effective zenith angle corresponding to the diurnally averaged insolation
and earth-sun distance at a particular latitude given the date.

`param_set` is an AbstractParameterSet from CLIMAParameters.jl.

`eot_correction` is an optional Boolean keyword argument that defaults to true
when set to true the equation of time correction is turned on.
This switch functionality is implemented for easy comparisons with reanalyses.

`milankovitch` is an optional Boolean keyword argument that defaults to true
when set to true the orbital parameters are calculated for the given DateTime,
when set to false the orbital parameters at the J2000 epoch from CLIMAParameters are used.
"""
function daily_zenith_angle(date::DateTime,
                            od::OrbitalData,
                            latitude::FT,
                            param_set::IP.AIP;
                            eot_correction::Bool=true,
                            milankovitch::Bool=true) where {FT}
    ϕ = deg2rad(latitude)

    epoch_string = "2000-01-01T11:58:56.816"
    # epoch_string::String = IP.epoch(param_set)
    Δt_years = get_Δt_years(param_set, date, epoch_string)
    
    # calculate orbital parameters or take values at J2000
    (ϖ, γ, e) = milankovitch ? compute_orbital_parameters(od, Δt_years) : compute_orbital_parameters(param_set)

    d, δ, _ = distance_declination_hourangle(date, (ϖ, γ, e), epoch_string, param_set, eot_correction)

    # sunrise/sunset angle
    T = tan(ϕ) * tan(δ)
    if T >= FT(1)
        ηd = FT(π)
    elseif T <= FT(-1)
        ηd = FT(0)
    else
        ηd = acos(FT(-1) * T)
    end

    # effective zenith angle to get diurnally averaged insolation (i.e., averaging cosine of zenith angle)
    daily_θ = mod(
        acos(FT(1 / π) * (ηd * sin(ϕ) * sin(δ) + cos(ϕ) * cos(δ) * sin(ηd))),
        FT(2π),
    )

    return daily_θ, d
end

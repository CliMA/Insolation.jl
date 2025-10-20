export instantaneous_zenith_angle, daily_zenith_angle

# true anomaly: angular distance from perihelion [radians], accurate to 
# O(e^4) where e is the eccentricity (see Fitzpatrick (2012), appendix A.10))
function true_anomaly(MA::FT, e::FT) where {FT <: Real}
    # Series expansion for true anomaly
    TA = MA + (2 * e - FT(1 / 4) * e^3) * sin(MA)
    + FT(5 / 4) * e^2 * sin(2 * MA)
    + FT(13 / 12) * e^3 * sin(3 * MA)
    return mod(TA, FT(2π))
end

# equation of time [radians]; this value can be scaled by `day_length / (2π)` 
# to get a time correction in seconds.
function equation_of_time(e::FT, MA::FT, γ::FT, ϖ::FT) where {FT <: Real}
    _Δt = -2 * e * sin(MA) + tan(γ / 2)^2 * sin(2 * (MA + ϖ))
    return mod(_Δt + FT(π), FT(2π)) - FT(π)
end

# calculate orbital parameters (Milankovitch)
compute_orbital_parameters(od, Δt_years::FT) where {FT} = (
    FT(od.ϖ_spline(Δt_years)),
    FT(od.γ_spline(Δt_years)),
    FT(od.e_spline(Δt_years)),
)

compute_orbital_parameters(param_set) = (
    IP.lon_perihelion_epoch(param_set),
    IP.obliq_epoch(param_set),
    IP.eccentricity_epoch(param_set),
)

function get_Δt_years(
    param_set::IP.InsolationParameters{FT},
    date::DateTime,
    time_of_epoch::DateTime,
) where {FT}
    (; year_anom, day) = param_set
    days_per_year = year_anom / day
    return FT(datetime2julian(date) - datetime2julian(time_of_epoch)) /
           days_per_year
end

"""
    _compute_distance_and_declination(
        Δt_years::FT,
        (ϖ, γ, e)::Tuple{FT, FT, FT},
        param_set::IP.AIP,
    ) where {FT}

Internal helper to compute Earth-Sun distance and declination angle.
This is called by both daily and instantaneous calculations.

Returns tuple `(d, δ, MA)`:
- `d`: Earth-Sun distance [m]
- `δ`: Declination angle [radians]
- `MA`: Mean anomaly [radians] (needed for hour angle)
"""
function _compute_distance_and_declination(
    Δt_years::FT,
    (ϖ, γ, e)::Tuple{FT, FT, FT},
    param_set::IP.AIP,
) where {FT}
    Ya = IP.year_anom(param_set)
    day_length = IP.day(param_set)
    d0 = IP.orbit_semimaj(param_set)
    M0 = IP.mean_anom_epoch(param_set)

    # mean anomaly given mean anomaly at epoch M0
    MA = mod(FT(2π) * (Δt_years) + M0, FT(2π))

    # true anomaly [radians]
    TA = true_anomaly(MA, e)

    # true longitude [radians]
    TL = mod(TA + ϖ, FT(2π))

    # declination [radians]
    δ = mod(asin(sin(γ) * sin(TL)), FT(2π))

    # earth-sun distance [m] 
    d = d0 * (1 - e^2) / (1 + e * cos(TA))

    return d, δ, MA
end

"""
    distance_declination_hourangle(
        Δt_years::FT,
        date::DateTime,
        time_of_epoch::DateTime,
        (ϖ, γ, e)::Tuple{FT, FT, FT},
        param_set::IP.AIP,
        eot_correction::Bool,
    ) where {FT}

Returns the Earth-Sun distance [m], declination angle [radians] and
hour angle [radians] at 0° longitude.

Used for instantaneous calculations.
"""
function distance_declination_hourangle(
    Δt_years::FT,
    date::DateTime,
    time_of_epoch::DateTime,
    (ϖ, γ, e)::Tuple{FT, FT, FT},
    param_set::IP.AIP,
    eot_correction::Bool,
) where {FT}
    Ya = IP.year_anom(param_set)
    day_length = IP.day(param_set)
    days_per_year = Ya / day_length
    
    d, δ, MA = _compute_distance_and_declination(Δt_years, (ϖ, γ, e), param_set)

    # hour angle, zero at local solar noon [radians]
    if eot_correction
        # radians to seconds
        Δt = equation_of_time(e, MA, γ, ϖ) / FT(2π) * day_length
    else
        Δt = FT(0)
    end
    time_in_day = FT(mod(datetime2julian(date), 1)) * day_length # s in this day
    η_UTC = mod(FT(2π) * (time_in_day + Δt) / day_length, FT(2π))

    return d, δ, η_UTC
end

"""
    instantaneous_zenith_angle(
        d::FT,
        δ::FT,
        η_UTC::FT,
        longitude::FT,
        latitude::FT,
    ) where {FT}

Returns the zenith angle, azimuthal angle, and Earth-Sun distance
at a particular longitude and latitude.

# Arguments
- `d::FT`: Earth-Sun distance [m]
- `δ::FT`: Declination angle [radians]
- `η_UTC::FT`: Hour angle at 0° longitude [radians]
- `longitude::FT`: Longitude [degrees]
- `latitude::FT`: Latitude [degrees]

# Returns
- `θ`: Solar zenith angle [radians]
- `ζ`: Solar azimuth angle [radians], 0 = due East, increasing CCW
- `d`: Earth-Sun distance [m]
"""
function instantaneous_zenith_angle(
    d::FT,
    δ::FT,
    η_UTC::FT,
    longitude::FT,
    latitude::FT,
) where {FT}
    λ = deg2rad(longitude)
    ϕ = deg2rad(latitude)
    # hour angle at longitude
    η = mod(η_UTC + λ, FT(2π))

    # zenith angle [radians]
    θ = mod(
        acos(
            max(FT(-1), min(FT(1), cos(ϕ) * cos(δ) * cos(η) + sin(ϕ) * sin(δ))),
        ),
        FT(2π),
    )

    # solar azimuth angle, ζ = 0 when due E and increasing CCW
    # ζ = 3π/2 (due S) when η=0 at local solar noon
    ζ = mod(
        FT(3π / 2) - atan(sin(η), cos(η) * sin(ϕ) - tan(δ) * cos(ϕ)),
        FT(2π),
    )

    # NOTE: d is returned here for API compatibility with InsolationCalc.jl,
    # even though it was passed in as an argument.
    return θ, ζ, d
end

function helper_instantaneous_zenith_angle(
    date::DateTime,
    od::OrbitalData,
    param_set::AIP;
    eot_correction = true,
)
    epoch_string = IP.epoch(param_set)
    time_of_epoch = DateTime(epoch_string, dateformat"y-m-dTHH:MM:SS.s")
    
    Δt_years = get_Δt_years(param_set, date, time_of_epoch)
    return distance_declination_hourangle(
        Δt_years, 
        date,
        time_of_epoch,
        Insolation.compute_orbital_parameters(od, Δt_years),
        param_set,
        eot_correction,
    )
end

function helper_instantaneous_zenith_angle(
    date::DateTime,
    param_set::AIP;
    eot_correction = true,
)
    epoch_string = IP.epoch(param_set)
    time_of_epoch = DateTime(epoch_string, dateformat"y-m-dTHH:MM:SS.s")

    Δt_years = get_Δt_years(param_set, date, time_of_epoch)
    return distance_declination_hourangle(
        Δt_years, 
        date,
        time_of_epoch,
        Insolation.compute_orbital_parameters(param_set),
        param_set,
        eot_correction,
    )
end

"""
    daily_zenith_angle(
        date::DateTime,
        od::OrbitalData,
        latitude::FT,
        param_set::IP.AIP;
        milankovitch::Bool = true
    ) where {FT <: Real}

Returns the effective zenith angle corresponding to the diurnally
averaged insolation, and the Earth-Sun distance.

# Arguments
- `date::DateTime`: Current date
- `od::OrbitalData`: Struct with orbital parameter splines
- `latitude::FT`: Latitude [degrees]
- `param_set::IP.AIP`: Parameter struct
- `milankovitch::Bool`: (default true) Use Milankovitch cycles. If false,
                         uses fixed epoch parameters.
"""
function daily_zenith_angle(
    date::DateTime,
    od::OrbitalData,
    latitude::FT,
    param_set::IP.AIP;
    milankovitch::Bool = true,
) where {FT}
    epoch_string = IP.epoch(param_set)
    time_of_epoch = DateTime(epoch_string, dateformat"y-m-dTHH:MM:SS.s")
    ϕ = deg2rad(latitude)

    Δt_years = get_Δt_years(param_set, date, time_of_epoch)

    # calculate orbital parameters or take values at J2000
    (ϖ, γ, e) =
        milankovitch ? compute_orbital_parameters(od, Δt_years) :
        compute_orbital_parameters(param_set)

    # Get distance and declination
    d, δ, _ = _compute_distance_and_declination(Δt_years, (ϖ, γ, e), param_set)

    # sunrise/sunset hour angle
    T = tan(ϕ) * tan(δ)
    if T >= FT(1)
        ηd = FT(π) # polar day
    elseif T <= FT(-1)
        ηd = FT(0) # polar night
    else
        ηd = acos(FT(-1) * T)
    end

    # effective zenith angle to get diurnally averaged insolation
    # (i.e., averaging cosine of zenith angle)
    daily_θ = mod(
        acos(FT(1 / π) * (ηd * sin(ϕ) * sin(δ) + cos(ϕ) * cos(δ) * sin(ηd))),
        FT(2π),
    )

    return daily_θ, d
end
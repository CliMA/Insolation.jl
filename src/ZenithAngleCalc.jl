export instantaneous_zenith_angle, daily_zenith_angle

# true anomaly, radians (eq 3.8)
function true_anomaly(MA::FT, e::FT) where {FT <: Real}
    TA = MA + (2*e - FT(1/4)*e^3)*sin(MA) 
            + FT(5/4)*e^2*sin(2*MA) 
            + FT(13/12)*e^3*sin(3*MA)
    return mod(TA, FT(2π))
end

# equation of time, radians
function equation_of_time(e::FT, MA::FT, γ::FT, ϖ::FT) where {FT <: Real}
    _Δt = -2*e*sin(MA) + tan(γ/2)^2*sin(2*(MA + ϖ));
    return mod(_Δt+FT(π), FT(2π)) - FT(π)
end

# calculate the distance, declination, and hour angle (at lon=0)
function distance_declination_hourangle(::Type{FT},
                                        date::DateTime,
                                        param_set::APS,
                                        eot_correction::Bool,
                                        milankovitch::Bool) where {FT <: Real}
    Ya::FT = year_anom(param_set)
    day_length::FT = Planet.day(param_set)
    AU::FT = astro_unit()
    _epoch::FT = epoch(param_set)
    M0::FT = mean_anom_epoch(param_set)

    t0_years = _epoch / Ya; # s --> years
    t_years = FT(datetime2julian(date)) * (day_length / Ya); # DateTime --> years
    
    # mean anomaly given mean anomaly at epoch (3.6)
    MA = mod(FT(2π) * (t_years - t0_years) + M0, FT(2π))

    # calculate orbital parameters or take values at J2000
    if milankovitch
        ϖ = FT(ϖ_spline(t_years - t0_years));
        γ = FT(γ_spline(t_years - t0_years));
        e = FT(e_spline(t_years - t0_years));
    else
        ϖ = FT(lon_perihelion_epoch(param_set));
        γ = FT(obliq_epoch(param_set));
        e = FT(eccentricity_epoch(param_set));
    end

    # true anomaly, radians (3.8)
    TA = true_anomaly(MA, e)

    # true longitude, radians (3.9)
    TL = mod(TA + ϖ, FT(2π))

    # declination, radians (3.16)
    δ = mod(asin(sin(γ) * sin(TL)), FT(2π))

    # earth-sun distance, (3.1)
    d = AU * (1 - e^2) / (1 + e*cos(TA))

    # hour angle, zero at local solar noon, radians (3.17)
    if eot_correction
        Δt = equation_of_time(e, MA, γ, ϖ) / FT(2π) * day_length # radians to seconds
    else
        Δt = FT(0)
    end
    t0 = FT(mod(datetime2julian(date), 1.0)) * day_length; # s in this day
    η_UTC = mod(FT(2π) * (t0 + Δt) / day_length, FT(2π));

    return d, δ, η_UTC
end

"""
    instantaneous_zenith_angle(date::DateTime,
                               longitude::FT,
                               latitude::FT,
                               param_set::APS;
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
                                    longitude::FT,
                                    latitude::FT,
                                    param_set::APS; 
                                    eot_correction::Bool=true,
                                    milankovitch::Bool=true) where {FT <: Real}
    λ = deg2rad(longitude)
    ϕ = deg2rad(latitude)

    d, δ, η_UTC = distance_declination_hourangle(FT, date, param_set, eot_correction, milankovitch)

    # hour angle
    η = mod(η_UTC + λ, FT(2π))

    # zenith angle, radians (3.18)
    θ = mod(acos(cos(ϕ)*cos(δ)*cos(η) + sin(ϕ)*sin(δ)), FT(2π))

    # solar azimuth angle, ζ = 0 when due E and increasing CCW
    # ζ = 3π/2 (due S) when η=0 at local solar noon
    ζ = mod(FT(3π/2)- atan(sin(η), cos(η)*sin(ϕ) - tan(δ)*cos(ϕ)), FT(2π))

    return θ, ζ, d
end

"""
    daily_zenith_angle(date::DateTime,
                       latitude::FT,
                       param_set::APS;
                       eot_correction::Bool=true,
                       milankovitch::Bool=true) where {FT <: Real}

Returns the daily averaged zenith angle and earth-sun distance
at a particular latitude given the date and orbital parameters
obliquity, longitude of perihelion, and eccentricity
param_set is an AbstractParameterSet from CLIMAParameters.jl.

`eot_correction` is an optional Boolean keyword argument that defaults to true
when set to true the equation of time correction is turned on.
This switch functionality is implemented for easy comparisons with reanalyses.

`milankovitch` is an optional Boolean keyword argument that defaults to true
when set to true the orbital parameters are calculated for the given DateTime,
when set to false the orbital parameters at the J2000 epoch from CLIMAParameters are used.
"""
function daily_zenith_angle(date::DateTime,
                            latitude::FT,
                            param_set::APS;
                            eot_correction::Bool=true,
                            milankovitch::Bool=true) where {FT <: Real}
    ϕ = deg2rad(latitude)

    d, δ, _ = distance_declination_hourangle(FT, date, param_set, eot_correction,milankovitch)
    
    # sunrise/sunset angle (3.19)
    T = tan(ϕ) * tan(δ)
    if T >= FT(1)
        ηd = FT(π)
    elseif T <= FT(-1)
        ηd = FT(0)
    else
        ηd = acos(FT(-1)*T)
    end
    
    # daily averaged zenith angle (3.20)
    daily_θ = mod(acos(FT(1/π)*(ηd*sin(ϕ)*sin(δ) + cos(ϕ)*cos(δ)*sin(ηd))), FT(2π))

    return daily_θ, d
end

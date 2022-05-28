export instantaneous_zenith_angle, daily_zenith_angle
export orbital_params

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

"""
    orbital_params(dt::FT, param_set) where {FT <: Real}

This function returns the orbital parameters (ϖ, γ, e) at a given 
`dt` number of years after the `_epoch`.
The parameters vary due to Milankovitch cycles. 
The dominant 10 frequencies of these cycles are used based on a
Fourier analysis of the parameters as calculated in the
Lasker et al. (2004) paper.
Data from this paper are in the "src/data/INSOL.LA2004.BTL.csv" file.
# Berger A. and Loutre M.F. (1991) paper. 
# Data from this paper is in the "src/data/orbit91.tsv" file.
"""
function orbital_params(dt::FT) where {FT <: Real}
    x = np.loadtxt("INSOL.LA2004.BTL.csv", delimiter=",", skiprows=1)
    e = CubicSplineInterpolation(x[:,0], x[:,1])(dt)
    γ = CubicSplineInterpolation(x[:,0], x[:,2])(dt)
    ϖ = CubicSplineInterpolation(x[:,0], x[:,3])(dt)
    return [ϖ, γ, e]
end
# function orbital_params(dt::FT, param_set) where {FT <: Real}
#     ϖ0::FT = lon_perihelion_epoch(param_set)
#     γ0::FT = obliq_epoch(param_set)
#     e0::FT = eccentricity_epoch(param_set)

#     ϖ = FT(mod(ϖ0 + 2π*dt/(26e3) + 2π*dt/(112e3), 2π))
#     γ = FT(γ0 + deg2rad(1.2)*sin(2π*dt/(41e3)))
#     e = FT(e0 + 0.02*sin(2π*dt/(405e3)) + 0.01*sin(2π*dt/(110e3)))
#     return [ϖ, γ, e]
# end


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

    time = datetime2julian(date) * day_length # in seconds
    dt = FT(time - _epoch) / Ya
    
    # mean anomaly given mean anomaly at epoch (3.6)
    MA = mod(FT(2π) * dt + M0, FT(2π))

    # calculate orbital parameters or take values at J2000
    if milankovitch
        ϖ, γ, e = orbital_params(dt, param_set)
    else
        ϖ::FT = lon_perihelion_epoch(param_set)
        γ::FT = obliq_epoch(param_set)
        e::FT = eccentricity_epoch(param_set)
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
    η_UTC = mod(FT(2π) * FT(time + Δt) / day_length, FT(2π))

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
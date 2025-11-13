"""
    Insolation

A Julia package to calculate top-of-atmosphere (TOA) insolation
(incoming solar radiation) based on Earth's (or another planetary 
body's) orbital parameters.

The calculations follow fundamental principles of celestial mechanics
and solar geometry, as described in ["Physics of Earth's Climate"](https://climate-dynamics.org/wp-content/uploads/2017/04/Climate_Book.pdf) 
by Tapio Schneider and Lenka Novak.

The package provides functions to:
- Calculate instantaneous insolation for a specific time and location.
- Calculate diurnally averaged insolation.
- Fetch and use orbital parameters (eccentricity, obliquity, and
  longitude of perihelion) from Laskar et al. (2004) 
  to compute insolation for paleoclimate studies.
"""
module Insolation

using Dates, DelimitedFiles, Interpolations
using Artifacts
using Adapt

include("Parameters.jl")
import .Parameters as IP
import .Parameters: InsolationParameters
const AIP = IP.AbstractInsolationParams

export orbital_params,
    OrbitalDataSplines, insolation, daily_insolation, TSIDataSpline, evaluate
export InsolationParameters

"""
    OrbitalDataSplines

A container struct that holds cubic spline interpolators for Earth's
orbital parameters, based on the Laskar 2004 dataset.

The time series of orbital parameters are lazily downloaded from the 
`orbital_parameters_dataset_path(artifact_dir)` path where `artifact_dir` is 
the path and filename to save the artifacts toml file.

The splines are functions of time (in years since J2000 epoch).

# GPU Support
This struct is GPU-compatible via Adapt.jl. To transfer to GPU memory:
```julia
using CUDA, Adapt
cpu_od = OrbitalDataSplines()  # Create on CPU
gpu_od = adapt(CuArray, cpu_od)  # Transfer to GPU
```
"""
struct OrbitalDataSplines{E,G,O}
    "Spline for eccentricity (e) [unitless]"
    e_spline::E
    "Spline for obliquity (γ) [radians]"
    γ_spline::G
    "Spline for longitude of perihelion (ϖ) [radians]"
    ϖ_spline::O
end

# Make OrbitalDataSplines GPU-compatible via Adapt.jl
# This allows the splines to be transferred to GPU memory using:
# gpu_od = adapt(CuArray, cpu_od)
Adapt.@adapt_structure OrbitalDataSplines

function OrbitalDataSplines()
    datapath = joinpath(artifact"laskar2004", "INSOL.LA2004.BTL.csv")
    Tx = Tuple{Matrix{Float64},Matrix{AbstractString}}
    laskar_data, _ = readdlm(datapath, ',', Float64, header = true)::Tx

    # Create a time range in years, with a 1000-year (1 kyr) step
    t_range = ((laskar_data[1, 1]*1000):1000:(laskar_data[end, 1]*1000))

    e_spline =
        cubic_spline_interpolation(t_range, laskar_data[:, 2]; extrapolation_bc = NaN)
    γ_spline =
        cubic_spline_interpolation(t_range, laskar_data[:, 3]; extrapolation_bc = NaN)
    ϖ_spline =
        cubic_spline_interpolation(t_range, laskar_data[:, 4]; extrapolation_bc = NaN)

    E = typeof(e_spline)
    G = typeof(γ_spline)
    O = typeof(ϖ_spline)
    return OrbitalDataSplines{E,G,O}(e_spline, γ_spline, ϖ_spline)
end

Base.broadcastable(x::OrbitalDataSplines) = tuple(x)

"""
    orbital_params(od::OrbitalDataSplines, dt::FT) where {FT <: Real}

Interpolates time-varying orbital parameters using Milankovitch cycles.

Interpolates orbital parameters from the Laskar et al. (2004) dataset for
paleoclimate studies. The parameters vary over geological timescales.

# Arguments
- `od::OrbitalDataSplines`: The struct containing orbital parameter splines.
- `dt::FT`: The time for interpolation [Years since J2000 epoch].

# Returns
- `(ϖ, γ, e)`: A tuple containing:
  - `ϖ`: Longitude of perihelion [radians]
  - `γ`: Obliquity (axial tilt) [radians]
  - `e`: Orbital eccentricity [unitless]
"""
function orbital_params(od::OrbitalDataSplines, dt::FT) where {FT<:Real}
    # Call the spline fields directly
    ϖ = od.ϖ_spline(dt)
    γ = od.γ_spline(dt)
    e = od.e_spline(dt)
    return ϖ, γ, e
end

"""
    orbital_params(param_set::AIP)

Returns fixed orbital parameters at epoch.

Uses the constant orbital parameter values at the reference epoch (typically
J2000) stored in the parameter set. Suitable for modern climate simulations
where orbital variations are negligible.

# Arguments
- `param_set::AIP`: Parameter struct containing epoch orbital parameters

# Returns
- `(ϖ, γ, e)`: A tuple containing:
  - `ϖ`: Longitude of perihelion at epoch [radians]
  - `γ`: Obliquity (axial tilt) at epoch [radians]
  - `e`: Orbital eccentricity at epoch [unitless]
"""
function orbital_params(param_set::AIP)
    ϖ = IP.lon_perihelion_epoch(param_set)
    γ = IP.obliq_epoch(param_set)
    e = IP.eccentricity_epoch(param_set)
    return ϖ, γ, e
end

"""
    TSIDataSpline

A container struct that holds a cubic interpolator for the monthly mean total
solar irradiance.

The spline is a function of date between `1850-01-15T12:00` and
`2229-12-15T12:00`. Dates outside this range will result in `NaN`.

# GPU Support
This struct is GPU-compatible via Adapt.jl. To transfer to GPU memory:
```julia
using CUDA, Adapt
cpu_tsi = TSIDataSpline(Float32)  # Create on CPU
gpu_tsi = adapt(CuArray, TSIDataSpline)  # Transfer to GPU
```
"""
struct TSIDataSpline{SPLINE,DATES,DATE<:Dates.DateTime}
    tsi_spline::SPLINE
    dates::DATES
    ref_date::DATE
end

"""
    _get_tsi_data()

Return the monthly dates and TSI data from the `cmip_monthly_tsi` artifact.
"""
function _get_tsi_data()
    datapath = joinpath(artifact"cmip_monthly_tsi", "cmip_monthly_tsi.csv")
    monthly_tsi_data, _ = readdlm(datapath, ',', String, header = true)
    monthly_dates = Dates.DateTime.(monthly_tsi_data[:, 1])
    tsi_data = parse.(Float32, monthly_tsi_data[:, 2])
    return monthly_dates, tsi_data
end

"""
    TSIDataSpline(::Type{FT}) where {FT <: AbstractFloat}

Constructs a `TSIDataSpline` that interpolates monthly mean total solar
irradiance as a function of the date and time.
"""
function TSIDataSpline(::Type{FT}) where {FT<:AbstractFloat}
    monthly_dates, tsi_data = _get_tsi_data()
    # Because times are not equispaced, we use the number of months since the
    # reference date as the data are monthly averages
    t_range = 0:(length(monthly_dates)-1)
    tsi_spline = cubic_spline_interpolation(t_range, FT.(tsi_data); extrapolation_bc = NaN)
    ref_date = first(monthly_dates)
    return TSIDataSpline(tsi_spline, monthly_dates, ref_date)
end

"""
    evaluate(tsi::TSIDataSpline, date::Dates.DateTime)

Interpolate monthly mean total solar irradiance at `date` via a cubic
interpolation.

The spline is a function of date between `1850-01-15T12:00` and
`2229-12-15T12:00`. Dates outside this range will result in `NaN`.
"""
function evaluate(tsi::TSIDataSpline, date::Dates.DateTime)
    if last(tsi.dates) == date
        return tsi.tsi_spline(length(tsi.dates) - 1)
    end
    if date < first(tsi.dates) || date > last(tsi.dates)
        return NaN
    end

    # There are two cases:
    # 1. Day and time of month is on or after 12:00 of the 15th
    # 2. Day and time of month is before 12:00 of the 15th
    day_val = Dates.day(date)
    hour_val = Dates.hour(date)
    if day_val >= 15 && hour_val >= 12
        # Cover casea 1
        target_year = Dates.year(date)
        target_month = Dates.month(date)
    else
        # Cover case 2
        target_month = Dates.month(date) - 1
        target_year = ifelse(target_month == 0, Dates.year(date) - 1, Dates.year(date))
        target_month = ifelse(target_month == 0, 12, target_month)
    end
    ref_year, ref_month = year(tsi.ref_date), month(tsi.ref_date)

    # Calculate whole months since reference date
    months_since_ref = (target_year - ref_year) * 12 + (target_month - ref_month)

    # Note that months_since_ref starts with zero
    # months_since_ref[i] is the reference date plus i months, so
    # tsi.dates[months_since_ref + 1] <= date <= tsi.dates[months_since_ref + 2]
    elapsed_milliseconds = date - tsi.dates[months_since_ref+1]
    total_milliseconds = tsi.dates[months_since_ref+2] - tsi.dates[months_since_ref+1]
    fractional_milliseconds = elapsed_milliseconds / total_milliseconds
    return tsi.tsi_spline(months_since_ref + fractional_milliseconds)
end

Adapt.@adapt_structure TSIDataSpline

Base.broadcastable(x::TSIDataSpline) = tuple(x)

include("SolarGeometry.jl")
include("InsolationCalc.jl")

end

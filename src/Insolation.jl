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

export orbital_params, OrbitalDataSplines, insolation, daily_insolation
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

include("SolarGeometry.jl")
include("InsolationCalc.jl")

end

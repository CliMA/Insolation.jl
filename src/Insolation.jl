"""
    Insolation

A Julia package to calculate top-of-atmosphere (TOA) insolation
(incoming solar radiation) based on Earth's orbital parameters.

The calculations follow fundamental principles of celestial mechanics
and solar geometry, as described in standard climate physics texts.

The package provides functions to:
- Calculate instantaneous insolation for a specific time and location.
- Calculate diurnally averaged insolation.
- Fetch and use orbital parameters (eccentricity, obliquity, and
  [cite_start]longitude of perihelion) from Laskar et al. (2004) [cite: 771]
  to compute insolation for paleoclimate studies.
"""
module Insolation

using Dates, DelimitedFiles, Interpolations
using Artifacts

include("Parameters.jl")
import .Parameters as IP
const AIP = IP.AbstractInsolationParams

export orbital_params,
    OrbitalData,
    insolation,
    solar_flux_and_cos_sza,
    daily_zenith_angle

"""
    OrbitalData

A container struct that holds cubic spline interpolators for Earth's
orbital parameters, based on the Laskar 2004 dataset.

The orbital parameters are lazily downloaded from the 
`orbital_parameters_dataset_path(artifact_dir)` path where `artifact_dir` is 
the path and filename to save the artifacts toml file.

The splines are functions of time (in years since J2000 epoch).
"""
struct OrbitalData{E, G, O}
    "Spline for eccentricity (e) [unitless]"
    e_spline::E
    "Spline for obliquity (γ) [radians]"
    γ_spline::G
    "Spline for longitude of perihelion (ϖ) [radians]"
    ϖ_spline::O
end

function OrbitalData()
    datapath = joinpath(artifact"laskar2004", "INSOL.LA2004.BTL.csv")
    Tx = Tuple{Matrix{Float64}, Matrix{AbstractString}}
    x, _ = readdlm(datapath, ',', Float64, header = true)::Tx
    
    # Create a time range in years, with a 1000-year (1 kyr) step
    t_range = ((x[1, 1] * 1000):1000:(x[end, 1] * 1000))
    
    e_spline =
        cubic_spline_interpolation(t_range, x[:, 2]; extrapolation_bc = NaN)
    γ_spline =
        cubic_spline_interpolation(t_range, x[:, 3]; extrapolation_bc = NaN)
    ϖ_spline =
        cubic_spline_interpolation(t_range, x[:, 4]; extrapolation_bc = NaN)

    E = typeof(e_spline)
    G = typeof(γ_spline)
    O = typeof(ϖ_spline)
    return OrbitalData{E, G, O}(e_spline, γ_spline, ϖ_spline)
end

Base.broadcastable(x::OrbitalData) = tuple(x)

"""
    orbital_params(od::OrbitalData, dt::FT) where {FT <: Real}

Interpolates orbital parameters for a given time `dt` using the
splines in `od`.

# Arguments
- `od::OrbitalData`: The struct containing orbital parameter splines.
- `dt::FT`: The time for interpolation [Years since J2000 epoch].

# Returns
- `(ϖ, γ, e)`: A tuple containing the longitude of perihelion [radians],
  obliquity [radians], and eccentricity [unitless].
"""
function orbital_params(od::OrbitalData, dt::FT) where {FT <: Real}
    # Call the spline fields directly
    ϖ = od.ϖ_spline(dt)
    γ = od.γ_spline(dt)
    e = od.e_spline(dt)
    return ϖ, γ, e
end

include("ZenithAngleCalc.jl")
include("InsolationCalc.jl")

# For backwards compatibility with package extensions
if !isdefined(Base, :get_extension)
    include(joinpath("..", "ext", "CreateParametersExt.jl"))
end

end # module

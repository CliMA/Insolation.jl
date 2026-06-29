export insolation, daily_insolation

"""
    solar_flux(
        d::FT,
        param_set::IP.AIP,
        date::Union{DateTime,Nothing} = nothing,
        solar_variability_spline::Union{TSIDataSpline,Nothing} = nothing,
    ) where {FT<:Real}

Calculate the solar radiative energy flux at the top of the atmosphere
(TOA) based on the planet-star distance and the inverse square law.

The total solar irradiance `S0` (either `tot_solar_irrad` or, if a
`solar_variability_spline` is supplied, the interpolated value) is defined at the
*mean* orbital distance (the semi-major axis). Since `d` is the planet-star
distance in units of the semi-major axis, the flux is simply
``S = S_0 / d^2``. The absolute size of the orbit does not enter: the flux
depends only on `S0` and the orbital shape (eccentricity and true anomaly).

# Arguments
- `d::FT`: Planet-star distance, in units of the semi-major axis [unitless]
- `param_set::IP.AIP`: Struct containing `tot_solar_irrad`, the total solar irradiance
  at the mean orbital distance [W m⁻²]
- `date::Union{DateTime,Nothing}`: (default `nothing`) Current date and time used to evaluate
  the solar flux if `solar_variability_spline` is available.
- `solar_variability_spline::Union{TSIDataSpline,Nothing}`: (default `nothing`) Use time-varying
  solar luminosity if a `TSIDataSpline` is passed as an argument.

# Returns
- `S`: Solar flux at the given planet-star distance [W m⁻²]
"""
function solar_flux(
    d::FT,
    param_set::IP.AIP,
    date::Union{DateTime,Nothing} = nothing,
    solar_variability_spline::Union{TSIDataSpline,Nothing} = nothing,
) where {FT<:Real}
    S0::FT =
        isnothing(solar_variability_spline) ? IP.tot_solar_irrad(param_set) :
        FT(evaluate(solar_variability_spline, date))

    # Solar radiative energy flux. `d` is in units of the semi-major axis (the mean
    # orbital distance at which S0 is defined), so the inverse-square law is S0 / d².
    S = S0 / d^2
    return S
end

"""
    insolation(
        θ::FT,
        d::FT,
        param_set::IP.AIP,
        date::Union{DateTime,Nothing} = nothing,
        solar_variability_spline::Union{TSIDataSpline,Nothing} = nothing,
    ) where {FT <: Real}

Calculate top-of-atmosphere (TOA) insolation and the cosine of the solar zenith angle.

Implements ``F = S \\cos(\\theta)`` where S is the solar flux at the given
planet-star distance. Insolation is set to 0 at night (when ``\\cos(\\theta) < 0``).

# Arguments
- `θ::FT`: Solar zenith angle [radians]
- `d::FT`: Planet-star distance, in units of the semi-major axis [unitless]
- `param_set::IP.AIP`: Parameter struct
- `date::Union{DateTime,Nothing}`: (default `nothing`) Current date and time used to evaluate the
  solar flux if `solar_variability_spline` is available.
- `solar_variability_spline::Union{TSIDataSpline,Nothing}`: (default `nothing`) Use time-varying
  solar luminosity if a `TSIDataSpline` is passed as an argument.

# Returns
A `NamedTuple` with fields:
- `F`: TOA insolation [W m⁻²]
- `S`: Solar flux at the given planet-star distance [W m⁻²]
- `μ`: Cosine of solar zenith angle [unitless], clamped to [0, 1]
"""
function insolation(
    θ::FT,
    d::FT,
    param_set::IP.AIP,
    date::Union{DateTime,Nothing} = nothing,
    solar_variability_spline::Union{TSIDataSpline,Nothing} = nothing,
) where {FT<:Real}
    # Calculate solar radiative energy flux (W m⁻²)
    S = solar_flux(d, param_set, date, solar_variability_spline)

    # Cosine of solar zenith angle (set to 0 at night)
    μ = max(FT(0), cos(θ))

    # TOA insolation
    F = S * μ

    return (; F, S, μ)
end

"""
    insolation(
        date::DateTime,
        latitude::FT1,
        longitude::FT2,
        param_set::IP.AIP,
        orbital_data::Union{OrbitalDataSplines, Nothing} = nothing,
        milankovitch::Bool = false,
        solar_variability_spline::Union{TSIDataSpline, Nothing} = nothing,
        eot_correction::Bool = true,
    ) where {FT1 <: Real, FT2 <: Real}

Calculate instantaneous TOA insolation with optional long-term variations
in Earth's orbital parameters (Milankovitch cycles) and solar luminosity.

# Arguments
- `date::DateTime`: Current date and time
- `latitude::FT1`: Latitude [degrees]
- `longitude::FT2`: Longitude [degrees]
- `param_set::IP.AIP`: Parameter struct
- `orbital_data::Union{OrbitalDataSplines, Nothing}`: (default nothing) Orbital parameter splines.
  **Required** when `milankovitch=true` for GPU compatibility.
- `milankovitch::Bool`: (default false) Use time-varying orbital parameters (Milankovitch cycles).
  The `OrbitalDataSplines` are Earth-specific (Laskar et al., 2004), so this is meaningful for Earth only.
- `solar_variability_spline::Union{TSIDataSpline, Nothing}`: (default nothing) Use time-varying
  solar luminosity if `TSIDataSpline` is passed as an argument. The `TSIDataSpline` gives the
  Sun's irradiance at 1 au, so this is meaningful for Earth only.
- `eot_correction::Bool`: (default true) Apply equation of time correction

# Returns
A `NamedTuple` with fields:
- `F`: TOA insolation [W m⁻²]
- `S`: Solar flux [W m⁻²]
- `μ`: Cosine of solar zenith angle [unitless]
- `ζ`: Solar azimuth angle [radians], 0 = due East, increasing counterclockwise

# Examples
```julia
# Modern climate (fixed epoch parameters)
(; F, S, μ, ζ) = insolation(date, lat, lon, param_set)

# Paleoclimate with Milankovitch cycles
od = OrbitalDataSplines()  # Load once
milankovitch = true
(; F, S, μ, ζ) = insolation(date, lat, lon, param_set, od, milankovitch)

# Without equation of time correction
orbital_data = nothing
milankovitch = false
solar_variability_spline = nothing
eot_correction = false
result = insolation(date, lat, lon, param_set, orbital_data, milankovitch, solar_variability_spline, eot_correction)
```

See also [`daily_insolation`](@ref) and [`solar_geometry`](@ref).

# GPU Usage
For GPU execution, create orbital and solar variability data on CPU and transfer
to GPU using Adapt.jl:
```julia
using CUDA, Adapt
cpu_od = OrbitalDataSplines()  # Create on CPU
gpu_od = adapt(CuArray, cpu_od)  # Transfer to GPU
cpu_solar = TSIDataSpline(Float32) # Create on CPU
gpu_solar = adapt(CuArray, cpu_solar) # Transfer to GPU
# In GPU kernel:
milankovitch=true
result = insolation(date, lat, lon, param_set, gpu_od, milankovitch, gpu_solar)
```
"""
function insolation(
    date::DateTime,
    latitude::FT1,
    longitude::FT2,
    param_set::IP.AIP,
    orbital_data::Union{OrbitalDataSplines,Nothing} = nothing,
    milankovitch::Bool = false,
    solar_variability_spline::Union{TSIDataSpline,Nothing} = nothing,
    eot_correction::Bool = true,
) where {FT1<:Real,FT2<:Real}
    # Get orbital parameters using helper function
    orb_params = get_orbital_parameters(
        date,
        param_set,
        orbital_data,
        milankovitch,
        eltype(param_set),
    )

    # Get solar geometry
    (; d, θ, ζ) = Insolation.solar_geometry(
        date,
        latitude,
        longitude,
        orb_params,
        param_set;
        eot_correction,
    )

    # Calculate insolation
    (; F, S, μ) = insolation(θ, d, param_set, date, solar_variability_spline)

    return (; F, S, μ, ζ)
end

"""
    daily_insolation(
        date::DateTime,
        latitude::Real,
        param_set::IP.AIP,
        orbital_data::Union{OrbitalDataSplines, Nothing} = nothing,
        milankovitch::Bool = false,
        solar_variability_spline::Union{TSIDataSpline, Nothing} = nothing,
    )

Calculate diurnally averaged TOA insolation with optional long-term variations
in orbital parameters (Milankovitch cycles) and solar luminosity. The insolation is
averaged over a full day.

# Arguments
- `date::DateTime`: Current date
- `latitude::Real`: Latitude [degrees]
- `param_set::IP.AIP`: Parameter struct
- `orbital_data::Union{OrbitalDataSplines, Nothing}`: (default nothing) Orbital parameter splines.
  **Required** when `milankovitch=true` for GPU compatibility.
- `milankovitch::Bool`: (default false) Use time-varying orbital parameters (Milankovitch cycles).
  The `OrbitalDataSplines` are Earth-specific (Laskar et al., 2004), so this is meaningful for Earth only.
- `solar_variability_spline::Union{TSIDataSpline, Nothing}`: (default nothing) Use time-varying
  solar luminosity if `TSIDataSpline` is passed as an argument. The `TSIDataSpline` gives the
  Sun's irradiance at 1 au, so this is meaningful for Earth only.

# Returns
A `NamedTuple` with fields:
- `F`: Daily averaged TOA insolation [W m⁻²]
- `S`: Solar flux [W m⁻²]
- `μ`: Daily averaged cosine of solar zenith angle [unitless]

# Examples
```julia
# Modern climate (fixed epoch parameters)
result = daily_insolation(date, lat, param_set)
# Access fields: result.F, result.S, result.μ

# Paleoclimate with Milankovitch cycles
od = OrbitalDataSplines()  # Load once
milankovitch = true
result = daily_insolation(date, lat, param_set, od, milankovitch)
```

See also [`insolation`](@ref).

# GPU Usage
For GPU execution, create orbital and solar variability data on CPU and transfer
to GPU using Adapt.jl:
```julia
using CUDA, Adapt
cpu_od = OrbitalDataSplines()  # Create on CPU
gpu_od = adapt(CuArray, cpu_od)  # Transfer to GPU
cpu_solar = TSIDataSpline(Float32) # Create on CPU
gpu_solar = adapt(CuArray, cpu_solar)
milankovitch = true
# In GPU kernel:
result = daily_insolation(date, lat, param_set, gpu_od, milankovitch, gpu_solar)
```
"""
function daily_insolation(
    date::DateTime,
    latitude::Real,
    param_set::IP.AIP,
    orbital_data::Union{OrbitalDataSplines,Nothing} = nothing,
    milankovitch::Bool = false,
    solar_variability_spline::Union{TSIDataSpline,Nothing} = nothing,
)
    # Get orbital parameters using helper function
    orb_params = get_orbital_parameters(
        date,
        param_set,
        orbital_data,
        milankovitch,
        eltype(param_set),
    )

    # Get effective zenith angle and distance for daily averaged insolation
    # (daily_distance_zenith_angle converts inputs to eltype(param_set) internally)
    (; daily_θ, d) =
        Insolation.daily_distance_zenith_angle(date, latitude, orb_params, param_set)

    # Return daily averaged insolation
    return insolation(daily_θ, d, param_set, date, solar_variability_spline)
end

"""
    get_orbital_parameters(
        date::DateTime,
        param_set::IP.AIP,
        orbital_data::Union{OrbitalDataSplines, Nothing},
        milankovitch::Bool,
        ::Type{FT},
    ) where {FT <: AbstractFloat}

Helper function to get orbital parameters with optional Milankovitch cycles.

Returns a tuple (ϖ, γ, e) of orbital parameters, selecting between epoch values
and time-varying Milankovitch values based on the `milankovitch` flag.

# Arguments
- `date::DateTime`: Current date
- `param_set::IP.AIP`: Parameter struct
- `orbital_data::Union{OrbitalDataSplines, Nothing}`: Pre-loaded orbital data
- `milankovitch::Bool`: Whether to use time-varying parameters
- `FT::Type`: Floating-point type

# Returns
- `(ϖ, γ, e)::Tuple{FT, FT, FT}`: Orbital parameters
"""
function get_orbital_parameters(
    date::DateTime,
    param_set::IP.AIP,
    orbital_data::Union{OrbitalDataSplines,Nothing},
    milankovitch::Bool,
    ::Type{FT},
) where {FT<:AbstractFloat}
    # Compute time-varying parameters if needed
    if milankovitch
        # Require pre-loaded orbital data for GPU compatibility
        if isnothing(orbital_data)
            error(
                "Spline interpolator orbital_data must be provided when milankovitch=true for GPU compatibility.\n
                Load OrbitalDataSplines: od = OrbitalDataSplines();\n
                Transfer to GPU: gpu_od = adapt(CuArray, od);\n
                Then call: insolation(date, lat, lon, param_set, gpu_od, true)\n",
            )
        end
        # The Laskar (2004) tables are indexed in Julian years since J2000
        Δt_years = Insolation.julian_years_since_epoch(param_set, date)
        ϖ, γ, e = Insolation.orbital_params(orbital_data, Δt_years)
    else
        # Compute epoch parameters
        ϖ, γ, e = Insolation.orbital_params(param_set)
    end

    return (FT(ϖ), FT(γ), FT(e))
end

export insolation, daily_insolation

"""
    solar_flux(
        d::FT,
        param_set::IP.AIP,
        date::Union{DateTime,Nothing} = nothing,
        solar_variability_spline::Union{TSIDataSpline,Nothing} = nothing,
    ) where {FT<:Real}

Calculates the solar radiative energy flux at the top of the atmosphere
(TOA) based on the planet-star distance and the inverse square law.

# Arguments
- `d::FT`: Planet-star distance [m]
- `param_set::IP.AIP`: Struct containing `tot_solar_irrad` [W m⁻²] and `orbit_semimaj` [m]
- `date::Union{DateTime,Nothing}`: (default `nothing`) Current date and time to evaluate the
  solar flux if `tsi_spline` is available.
- `solar_variability_spline::Union{TSIDataSpline,Nothing}`: Use time-varying solar
  luminosity if `TSIDataSpline` is passed as an argument.
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
    d0::FT = IP.orbit_semimaj(param_set)

    # Solar radiative energy flux
    S = S0 * (d0 / d)^2
    return S
end

"""
    insolation(θ::FT, d::FT, param_set::IP.AIP) where {FT <: Real}

Calculates top-of-atmosphere (TOA) insolation and cosine of solar zenith angle.

Implements ``F = S \\cos(\\theta)`` where S is the solar flux at the given
planet-star distance. Insolation is set to 0 at night (when ``\\cos(\\theta) < 0``).

# Arguments
- `θ::FT`: Solar zenith angle [radians]
- `d::FT`: Planet-star distance [m]
- `param_set::IP.AIP`: Parameter struct
- `date::Union{DateTime,Nothing}`: (default `nothing`) Current date and time to evaluate the
  solar flux if `tsi_spline` is available.
- `solar_variability_spline::Union{TSIDataSpline,Nothing}`: Use time-varying solar luminosity if
  `TSIDataSpline` is passed as an argument.

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

Calculates instantaneous TOA insolation with optional long-term variations
in Earth's orbital parameters (Milankovitch cycles) and solar luminosity.

# Arguments
- `date::DateTime`: Current date and time
- `latitude::FT1`: Latitude [degrees]
- `longitude::FT2`: Longitude [degrees]
- `param_set::IP.AIP`: Parameter struct
- `orbital_data::Union{OrbitalDataSplines, Nothing}`: (default nothing) Orbital parameter splines.
  **Required** when `milankovitch=true` for GPU compatibility.
- `milankovitch::Bool`: (default false) Use time-varying orbital parameters (Milankovitch cycles)
- `solar_variability_spline::Union{TSIDataSpline, Nothing}`: (default nothing) Use time-varying
  solar luminosity if `TSIDataSpline` is passed as an argument.
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
milankovitch = false,
solar_variability_spline = nothing
eot_correction = false
result = insolation(date, lat, lon, param_set, milankovitch, solar_variability_spline, eot_correction)
```

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

Calculates diurnally averaged TOA insolation with optional long-term variations
in orbital parameters (Milankovitch cycles) and solar luminosity. The insolation is 
averaged over a full day.

# Arguments
- `date::DateTime`: Current date
- `latitude::FT`: Latitude [degrees]
- `param_set::IP.AIP`: Parameter struct
- `orbital_data::Union{OrbitalDataSplines, Nothing}`: (default nothing) Orbital parameter splines.
  **Required** when `milankovitch=true` for GPU compatibility.
- `milankovitch::Bool`: (default false) Use time-varying orbital parameters (Milankovitch cycles).
- `solar_variability_spline::Union{TSIDataSpline, Nothing}`: (default nothing) Use time-varying
  solar luminosity if `TSIDataSpline` is passed as an argument.

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
result = daily_insolation(date, lat, param_set, od; milankovitch=true)
```

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
    (; daily_θ, d) = Insolation.daily_distance_zenith_angle(
        date,
        eltype(param_set)(latitude),
        orb_params,
        param_set,
    )

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
                Then call: insolation(date, lat, lon, param_set, gpu_od; milankovitch=true)\n",
            )
        end
        Δt_years = Insolation.years_since_epoch(param_set, date)
        ϖ, γ, e = Insolation.orbital_params(orbital_data, Δt_years)
    else
        # Compute epoch parameters
        ϖ, γ, e = Insolation.orbital_params(param_set)
    end

    return (FT(ϖ), FT(γ), FT(e))
end

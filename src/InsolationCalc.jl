export insolation, daily_insolation

"""
    solar_flux(d::FT, param_set::IP.AIP) where {FT <: Real}

Calculates the solar radiative energy flux at the top of the atmosphere
(TOA) based on the planet-star distance and the inverse square law.

# Arguments
- `d::FT`: Planet-star distance [m]
- `param_set::IP.AIP`: Struct containing `tot_solar_irrad` [W m⁻²] and `orbit_semimaj` [m]
"""
function solar_flux(d::FT, param_set::IP.AIP) where {FT <: Real}
    S0::FT = IP.tot_solar_irrad(param_set)
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

# Returns
- `F`: TOA insolation [W m⁻²]
- `S`: Solar flux at the given planet-star distance [W m⁻²]
- `μ`: Cosine of solar zenith angle [unitless], clamped to [0, 1]
"""
function insolation(θ::FT, d::FT, param_set::IP.AIP) where {FT <: Real}
    # Calculate solar radiative energy flux (W m⁻²)
    S = solar_flux(d, param_set)

    # Cosine of solar zenith angle (set to 0 at night) 
    μ = max(FT(0), cos(θ))

    # TOA insolation
    F = S * μ

    return F, S, μ
end

"""
    insolation(
        date::DateTime,
        latitude::FT1,
        longitude::FT2,
        param_set::IP.AIP,
        orbital_data::Union{OrbitalDataSplines, Nothing} = nothing;
        milankovitch::Bool = false,
        solar_variability::Bool = false,
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
- `solar_variability::Bool`: (default false) Use time-varying solar luminosity
- `eot_correction::Bool`: (default true) Apply equation of time correction

# Returns
- `F`: TOA insolation [W m⁻²]
- `S`: Solar flux [W m⁻²]
- `μ`: Cosine of solar zenith angle [unitless]
- `ζ`: Solar azimuth angle [radians], 0 = due East, increasing counterclockwise

# Examples
```julia
# Modern climate (fixed epoch parameters)
F, S, μ, ζ = insolation(date, lat, lon, param_set)

# Paleoclimate with Milankovitch cycles 
od = OrbitalDataSplines()  # Load once
F, S, μ, ζ = insolation(date, lat, lon, param_set, od; milankovitch=true)

# Without equation of time correction
F, S, μ, ζ = insolation(date, lat, lon, param_set; eot_correction=false)
```

# GPU Usage
For GPU execution, create orbital data on CPU and transfer to GPU using Adapt.jl:
```julia
using CUDA, Adapt
cpu_od = OrbitalDataSplines()  # Create on CPU
gpu_od = adapt(CuArray, cpu_od)  # Transfer to GPU
# In GPU kernel:
F, S, μ, ζ = insolation(date, lat, lon, param_set, gpu_od; milankovitch=true)
```
"""
function insolation(
    date::DateTime,
    latitude::FT1,
    longitude::FT2,
    param_set::IP.AIP,
    orbital_data::Union{OrbitalDataSplines, Nothing} = nothing;
    milankovitch::Bool = false,
    solar_variability::Bool = false,
    eot_correction::Bool = true,
) where {FT1 <: Real, FT2 <: Real}
    # Get orbital parameters using helper function
    orb_params = get_orbital_parameters(
        date,
        param_set,
        orbital_data,
        milankovitch,
        eltype(param_set),
    )

    # Get solar geometry
    d, θ, ζ = Insolation.solar_geometry(
        date,
        latitude,
        longitude,
        orb_params,
        param_set;
        eot_correction,
    )

    # Calculate insolation
    # Note: solar_variability is a placeholder for future solar luminosity variations
    F, S, μ = insolation(θ, d, param_set)

    return F, S, μ, ζ
end

"""
    daily_insolation(
        date::DateTime,
        latitude::Real,
        param_set::IP.AIP,
        orbital_data::Union{OrbitalDataSplines, Nothing} = nothing;
        milankovitch::Bool = false,
        solar_variability::Bool = false,
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
- `milankovitch::Bool`: (default false) Use time-varying orbital parameters (Milankovitch cycles)
- `solar_variability::Bool`: (default false) Use time-varying solar luminosity (placeholder for future)

# Returns
- `F`: Daily averaged TOA insolation [W m⁻²]
- `S`: Solar flux [W m⁻²]
- `μ`: Daily averaged cosine of solar zenith angle [unitless]

# Examples
```julia
# Modern climate (fixed epoch parameters)
F, S, μ = daily_insolation(date, lat, param_set)

# Paleoclimate with Milankovitch cycles
od = OrbitalDataSplines()  # Load once 
F, S, μ = daily_insolation(date, lat, param_set, od; milankovitch=true)
```

# GPU Usage
For GPU execution, create orbital data on CPU and transfer to GPU using Adapt.jl:
```julia
using CUDA, Adapt
cpu_od = OrbitalDataSplines()  # Create on CPU
gpu_od = adapt(CuArray, cpu_od)  # Transfer to GPU
# In GPU kernel:
F, S, μ = daily_insolation(date, lat, param_set, gpu_od; milankovitch=true)
```
"""
function daily_insolation(
    date::DateTime,
    latitude::Real,
    param_set::IP.AIP,
    orbital_data::Union{OrbitalDataSplines, Nothing} = nothing;
    milankovitch::Bool = false,
    solar_variability::Bool = false,
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
    daily_θ, d = Insolation.daily_distance_zenith_angle(
        date,
        eltype(param_set)(latitude),
        orb_params,
        param_set,
    )

    # Calculate daily averaged insolation
    # Note: solar_variability is a placeholder for future solar luminosity variations
    F, S, μ = insolation(daily_θ, d, param_set)

    return F, S, μ
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
    orbital_data::Union{OrbitalDataSplines, Nothing},
    milankovitch::Bool,
    ::Type{FT},
) where {FT <: AbstractFloat}
    # Compute time-varying parameters if needed
    if milankovitch
        # Require pre-loaded orbital data for GPU compatibility
        if isnothing(orbital_data)
            error(
                "Spline interpolator orbital_data must be provided when milankovitch=true for GPU compatibility. " *
                "Load OrbitalDataSplines: od = OrbitalDataSplines(); " *
                "Transfer to GPU: gpu_od = adapt(CuArray, od); " *
                "Then call: insolation(date, lat, lon, param_set, gpu_od; milankovitch=true)",
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

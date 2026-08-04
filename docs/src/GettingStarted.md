# Getting Started

This guide will help you get started with `Insolation.jl`, covering installation, basic usage, and common workflows.

## Installation

`Insolation.jl` is a registered Julia package. Install it using Julia's package manager:

```julia
using Pkg
Pkg.add("Insolation")
```

For the latest development version from GitHub:

```julia
Pkg.add(url="https://github.com/CliMA/Insolation.jl")
```

## Loading the Package

```julia
using Insolation
using Dates  # For DateTime objects
```

## Basic Concepts

### Coordinate System

- **Latitude**: Degrees North (positive) or South (negative), range [-90, 90]
- **Longitude**: Degrees East (positive) or West (negative), range [-180, 180]
- **Time**: Julia `DateTime` objects (UTC)

### Parameters

Insolation calculations require a parameter set containing orbital and physical constants. Create one using:

```julia
# Default Earth parameters (requires ClimaParams.jl)
using ClimaParams
params = InsolationParameters(Float64)

# Without ClimaParams, specify all parameters manually (the values below are
# the ClimaParams defaults for modern Earth):
params = InsolationParameters{Float64}(
    year_anom = 31558433.5,             # Anomalistic year [s] (≈ 365.2596 days)
    day = 86400.0,                      # Day length [s]
    eccentricity_epoch = 0.01671123,    # Eccentricity [unitless]
    obliq_epoch = deg2rad(23.43927944), # Obliquity [rad]
    lon_perihelion_epoch = deg2rad(282.93768193),  # Longitude of perihelion [rad]
    tot_solar_irrad = 1362.0,           # Solar irradiance at mean orbital distance [W m⁻²]
    epoch = DateTime(2000, 1, 1, 11, 58, 55, 816),  # J2000 epoch (UTC)
    mean_anom_epoch = deg2rad(357.52688973)  # Mean anomaly at epoch [rad]
)
```

The `epoch` is J2000, 12:00 Terrestrial Time on January 1, 2000, expressed in UTC. The
orbital elements are their values at that instant: the obliquity is the IAU 2006 mean
obliquity of the ecliptic, and the eccentricity, longitude of perihelion, and mean anomaly
are one self-consistent set from the JPL (Standish) 1800–2050 Earth-Moon barycenter
elements. Passing `milankovitch = true` replaces the three orbital elements with
time-varying values (see [Milankovitch Cycles](@ref)) but leaves `epoch`, `year_anom`,
`day`, and `mean_anom_epoch` untouched.

## Basic Usage

### Instantaneous Insolation

Calculate insolation at a specific time and location:

```@example howto
using Insolation
using Dates
using ClimaParams

# Setup
params = InsolationParameters(Float64)
date = DateTime(2024, 6, 21, 18, 0, 0)  # Summer solstice, 18h UTC
lat = 40.0    # Boulder, Colorado (degrees North)
lon = -105.0  # (degrees East)

# Calculate insolation and solar geometry
(; F, S, μ, ζ) = insolation(date, lat, lon, params)

println("TOA Insolation: $(F) W/m²")
println("Solar flux: $(S) W/m²")
println("Cosine of zenith angle: $(μ)")
println("Solar zenith angle: $(rad2deg(acos(μ)))°")
println("Solar azimuth angle: $(rad2deg(ζ))°")
```

### Daily-Averaged Insolation

Calculate diurnally averaged insolation (useful for climate models):

```@example howto
# Daily average only depends on date and latitude
date = DateTime(2024, 6, 21)
lat = 40.0

(; F, μ) = daily_insolation(date, lat, params)

println("Daily-averaged TOA Insolation: $(F) W/m²")
println("Daily-averaged cosine of zenith angle: $(μ)")
```

### Computing Solar Position

For applications that only need solar geometry without insolation:

```@example howto
# Get orbital parameters
orb_params = orbital_params(params)

# Calculate solar geometry
(; d, θ, ζ) = solar_geometry(date, lat, lon, orb_params, params)

println("Sun-Earth distance: $(d) (in units of the semi-major axis; ≈ AU for Earth)")
println("Solar zenith angle: $(rad2deg(θ))°")
println("Solar azimuth angle: $(rad2deg(ζ))°")
```

## Working with Milankovitch Cycles

For paleoclimate applications, use time-varying orbital parameters:

```@example howto
# Load orbital parameter time series (covers -50 to +20 Myr around J2000).
# Construction is expensive; do it once and reuse.
orbital_data = OrbitalDataSplines()

# Turning on Milankovitch cycles takes the orbital elements from the Laskar et al.
# (2004) tables at the requested date instead of from the parameter set
milankovitch = true
date = DateTime(-18000, 6, 21)  # ~20 kyr before J2000, in astronomical year numbering
F = insolation(date, lat, lon, params, orbital_data, milankovitch).F
println("TOA insolation 20 kyr ago: $(round(F, digits = 1)) W/m²")

# The orbital parameters themselves, indexed in Julian years since J2000
ϖ, γ, e = orbital_params(orbital_data, -20000.0)
println("Obliquity 20 kyr ago: $(round(rad2deg(γ), digits = 2))°")
println("Eccentricity 20 kyr ago: $(round(e, digits = 4))")
```

!!! warning "Calendar dates drift relative to the seasons"
    The date passed here fixes the position relative to *perihelion*, not relative to the
    equinoxes, so a given calendar date corresponds to a different season in deep time.
    See [Milankovitch Cycles](@ref) for the size of the drift and for
    calendar-independent ways to make paleoclimate comparisons.

## Advanced Options

### GPU Usage

`Insolation.jl` supports GPU computation through `Adapt.jl`. This is useful for large-scale climate simulations.

```julia
using CUDA, Adapt, Insolation, ClimaParams

# Create parameters and orbital data on CPU
params = InsolationParameters(Float32)
cpu_od = OrbitalDataSplines()
cpu_tsi = TSIDataSpline(Float32)

# Transfer orbital data to GPU
gpu_od = adapt(CuArray, cpu_od)
gpu_tsi = adapt(CuArray, cpu_tsi)

# Use in GPU kernels (positional argument required for GPU compatibility)
milankovitch = true
result = insolation(date, lat, lon, params, gpu_od, milankovitch, gpu_tsi)
```

**Design Note**: The constructor creates data on CPU; users explicitly transfer to GPU using `adapt()`. This gives explicit control over data placement and follows Julia GPU ecosystem conventions.

### Equation of Time Correction

Control whether to apply the equation of time correction:

```@example howto
# Without equation of time correction (simpler, less accurate; for testing against codes not using EoT corrections)
date = DateTime(2000, 1, 1, 13, 0, 0)
lat = 40.0
lon = 15.0

orbital_data = nothing
milankovitch = false
solar_variability_spline = nothing
eot_correction = false
result = insolation(date, lat, lon, params, orbital_data, milankovitch, solar_variability_spline, eot_correction)
```

### Custom Orbital Parameters

Override specific parameters for sensitivity studies:

```@example howto
# Obliquity raised 5° above the default of 23.44°
params_high_obliq = InsolationParameters(Float64, (;
    obliq_epoch = deg2rad(28.44)
))

F_modified = insolation(date, lat, lon, params_high_obliq).F
```

### Type Flexibility

Use Float32 for reduced memory usage (useful for large-scale simulations):

```@example howto
params_f32 = InsolationParameters(Float32)
lat_f32 = Float32(lat)
lon_f32 = Float32(lon)

result = insolation(date, lat_f32, lon_f32, params_f32)
```

## Common Patterns

### Loop Over Time and Space

```julia
using Insolation
using Dates
using ClimaParams

params = InsolationParameters(Float64)

# Create arrays
latitudes = -90.0:10.0:90.0
longitudes = -180.0:15.0:180.0
dates = DateTime(2024, 1, 1):Day(1):DateTime(2024, 12, 31)

# Allocate output
insolation_array = zeros(length(latitudes), length(longitudes), length(dates))

# Calculate
for (k, date) in enumerate(dates)
    for (j, lon) in enumerate(longitudes)
        for (i, lat) in enumerate(latitudes)
            F = insolation(date, lat, lon, params).F
            insolation_array[i, j, k] = F
        end
    end
end
```

### Seasonal Cycle at a Location

```julia
using Insolation
using Dates
using ClimaParams
using Plots

params = InsolationParameters(Float64)
lat = 40.0
lon = -105.0

# Sample every day for a year
dates = DateTime(2024, 1, 1):Day(1):DateTime(2024, 12, 31)
daily_insol = [daily_insolation(d, lat, params).F for d in dates]

# Plot
plot(dates, daily_insol,
     xlabel="Date",
     ylabel="Daily-mean TOA Insolation [W/m²]",
     title="Seasonal Cycle at $(lat)°N, $(lon)°E",
     legend=false)
```

## Next Steps

- See [Insolation Examples](@ref) for visualization and more complex use cases
- Learn about [Milankovitch Cycles](@ref) for paleoclimate applications
- Check the [API Reference](@ref) for complete function documentation

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

# Without ClimaParams, specify all parameters manually:
params = InsolationParameters{Float64}(
    year_anom = 365.259636 * 86400.0,  # Anomalistic year [s]
    day = 86400.0,                      # Day length [s]
    orbit_semimaj = 1.496e11,           # Semi-major axis [m]
    eccentricity_epoch = 0.0167,        # Eccentricity
    obliq_epoch = deg2rad(23.44),       # Obliquity [rad]
    lon_perihelion_epoch = deg2rad(282.95),  # Longitude of perihelion [rad]
    tot_solar_irrad = 1362.0,           # Solar irradiance at 1 AU [W/m²]
    epoch = DateTime(2000, 1, 1, 11, 58, 56, 816),  # J2000 epoch
    mean_anom_epoch = deg2rad(357.5291)  # Mean anomaly at epoch [rad]
)
```

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
F, S, μ, ζ = insolation(date, lat, lon, params)

println("TOA Insolation: $F W/m²")
println("Solar flux: $S W/m²")
println("Cosine of zenith angle: $μ")
println("Solar zenith angle: $(rad2deg(acos(μ)))°")
println("Solar azimuth angle: $(rad2deg(ζ))°")
```

### Daily-Averaged Insolation

Calculate diurnally averaged insolation (useful for climate models):

```@example howto
# Daily average only depends on date and latitude
date = DateTime(2024, 6, 21)
lat = 40.0

F_daily, S_daily, μ_daily = daily_insolation(date, lat, params)

println("Daily-averaged TOA Insolation: $F_daily W/m²")
println("Daily-averaged cosine of zenith angle: $μ_daily")
```

### Computing Solar Position

For applications that only need solar geometry without insolation:

```@example howto
# Get orbital parameters
orb_params = orbital_params(params)

# Calculate solar geometry
distance, zenith_angle, azimuth = solar_geometry(
    date, lat, lon, orb_params, params
)

println("Sun-Earth distance: $(distance / 1.496e11) AU")
println("Solar zenith angle: $(rad2deg(zenith_angle))°")
println("Solar azimuth angle: $(rad2deg(azimuth))°")
```

## Working with Milankovitch Cycles

For paleoclimate applications, use time-varying orbital parameters:

```@example howto
# Load orbital parameter time series (covers -50 to +20 million years around present)
orbital_data = OrbitalDataSplines()

# Calculate for 20,000 years ago (Last Glacial Maximum)
date = DateTime(2000, 1, 1)  # Reference date (day of year matters)
params = InsolationParameters(Float64)

# Use Milankovitch cycles
milankovitch = true
F, S, μ, ζ = insolation(
    date, lat, lon, params,
    orbital_data,
    milankovitch,
)

# Get orbital parameters at specific time
years_since_epoch = -20000.0  # 20 kyr ago
ϖ, γ, e = orbital_params(orbital_data, years_since_epoch)
println("Obliquity 20 kyr ago: $(rad2deg(γ))°")
println("Eccentricity 20 kyr ago: $e")
```

## Advanced Options

### GPU Usage

`Insolation.jl` supports GPU computation through `Adapt.jl`. This is useful for large-scale climate simulations.

```julia
using CUDA, Adapt, Insolation, ClimaParams

# Create parameters and orbital data on CPU
params = InsolationParameters(Float32)
cpu_od = OrbitalDataSplines()

# Transfer orbital data to GPU
gpu_od = adapt(CuArray, cpu_od)

# Use in GPU kernels (positional argument required for GPU compatibility)
milankovitch = true
F, S, μ, ζ = insolation(date, lat, lon, params, gpu_od, milankovitch)
```

**Design Note**: The constructor creates data on CPU; users explicitly transfer to GPU using `adapt()`. This gives explicit control over data placement and follows Julia GPU ecosystem conventions.

### Equation of Time Correction

Control whether to apply the equation of time correction:

```@example howto
# Without equation of time correction (simpler, less accurate; for testing against codes not using EoT corrections)
date = DateTime(2000, 1, 1, 13, 0, 0)
lat = 40.0
lon = 15.0

milankovitch = false
solar_variability = false
eot_correction = false
F, S, μ, ζ = insolation(date, lat, lon, params, milankovitch, solar_variability, eot_correction)
```

### Custom Orbital Parameters

Override specific parameters for sensitivity studies:

```@example howto
# Increase obliquity by 5 degrees
params_high_obliq = InsolationParameters(Float64, (;
    obliq_epoch = deg2rad(23.44 + 5.0)
))

F_modified, _, _, _ = insolation(date, lat, lon, params_high_obliq)
```

### Type Flexibility

Use Float32 for reduced memory usage (useful for large-scale simulations):

```@example howto
params_f32 = InsolationParameters(Float32)
lat_f32 = Float32(lat)
lon_f32 = Float32(lon)

F, S, μ, ζ = insolation(date, lat_f32, lon_f32, params_f32)
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
            F, _, _, _ = insolation(date, lat, lon, params)
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
daily_insol = [daily_insolation(d, lat, params)[1] for d in dates]

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

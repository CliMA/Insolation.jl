# Insolation.jl

```@meta
CurrentModule = Insolation
```

`Insolation.jl` is a Julia package for calculating solar radiation and solar geometry for climate science and Earth system modeling applications. It provides efficient and accurate computations of solar zenith angle, azimuth angle, and incoming solar radiation (insolation) at the top of atmosphere.

## Overview

The package computes:

- **Solar Geometry**: Zenith angle, azimuth angle, and planet-star distance for any location and time
- **Insolation**: Instantaneous and daily-averaged incoming solar radiation
- **Orbital Parameters**: Time-varying orbital elements using the [Laskar et al. (2004)](https://doi.org/10.1051/0004-6361:20041335) solution for paleoclimate applications

### Key Features

- **Planetary Flexibility**: Works for Earth by default, but applicable to any planetary body with configurable parameters
- **Milankovitch Cycles**: Includes time-varying orbital parameters for paleoclimate studies spanning from -50 to +20 Myr
- **Type-Generic**: Supports different floating-point types (Float32, Float64) for performance optimization
- **GPU-Compatible**: Designed with GPU execution in mind, avoiding dynamic allocations in computational kernels
- **CliMA Integration**: Built for use with [ClimaAtmos.jl](https://github.com/CliMA/ClimaAtmos.jl), [RRTMGP.jl](https://github.com/CliMA/RRTMGP.jl), and other CliMA packages

## Package Structure

The package is organized into several modules:

- **`Parameters.jl`**: Defines `InsolationParameters` struct containing physical and orbital parameters
- **`SolarGeometry.jl`**: Calculates solar geometry (zenith angle, azimuth, declination, distance)
- **`InsolationCalc.jl`**: Computes top-of-atmosphere insolation from solar geometry
- **`CreateParametersExt.jl`**: Extension for integration with [ClimaParams.jl](https://github.com/CliMA/ClimaParams.jl)

## Mathematical Background

The calculations follow fundamental principles of celestial mechanics and solar geometry, as described in [Physics of Earth's Climate](https://climate-dynamics.org/wp-content/uploads/2017/04/Climate_Book.pdf) by Tapio Schneider and Lenka Novak. See the [Mathematical Background](SolarGeometry.md) page for detailed mathematical formulations of the zenith angle, azimuth angle, and other astronomical calculations.

## Quick Example

```@example quick
using Pkg
Pkg.instantiate()  # Install all dependencies
using Insolation   
using Dates

# Calculate insolation for a specific location and time
FT = Float64
lat = FT(40.0)    # degrees North
lon = FT(-105.0)  # degrees East
date = DateTime(2024, 6, 21, 12, 0, 0)  # Summer solstice, noon

# Create Earth parameters
using ClimaParams
params = InsolationParameters(FT)

# Calculate instantaneous insolation with solar geometry
(; F, S, μ, ζ) = insolation(date, lat, lon, params)
# F: TOA insolation [W m⁻²]
# S: Solar flux [W m⁻²]  
# μ: Cosine of solar zenith angle
# ζ: Solar azimuth angle [radians]
```

This is the instantaneous insolation at the given location and time. Here's the daily averaged insolation at the same location and time:

```@example quick
# Calculate daily-averaged insolation
(; F, S, μ) = daily_insolation(date, lat, params)
```

## Documentation Outline

```@contents
Pages = [
    "GettingStarted.md",
    "InsolationExamples.md",
    "Milankovitch.md",
    "SolarGeometry.md",
    "library.md"
]
Depth = 1
```

## Authors

`Insolation.jl` is being developed by the [Climate Modeling Alliance](https://clima.caltech.edu). The package was created for use with [ClimaAtmos.jl](https://github.com/CliMA/ClimaAtmos.jl), [RRTMGP.jl](https://github.com/CliMA/RRTMGP.jl), and other CliMA packages as part of the broader CliMA ecosystem.

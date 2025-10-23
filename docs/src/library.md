# API Reference

```@docs
Insolation.Insolation
```

Complete documentation of the public API for `Insolation.jl`.

## Overview

The package provides functions organized into three main categories:

1. **Parameter Management**: Create and manage physical/orbital parameters
2. **Solar Geometry**: Calculate solar geometry (zenith angle, azimuth angle, distance)
3. **Insolation Calculations**: Compute top-of-atmosphere solar radiation

## Parameters

### Parameter Structures

**`InsolationParameters{FT}`** - The main parameter struct containing physical and orbital constants needed for insolation calculations. See the source code or API documentation for all fields.

### Parameter Creation

The `CreateParametersExt` extension provides convenient constructors when `ClimaParams.jl` is loaded. See the Extensions section for details.

When `ClimaParams.jl` is available, you can create parameters using:

```julia
# With ClimaParams.jl loaded
using Insolation
params = InsolationParameters(Float64)

# With custom overrides
params = InsolationParameters(Float64, (; tot_solar_irrad = 1365.0))
```

Without ClimaParams.jl, create parameters directly:

```julia
using Insolation.Parameters
params = InsolationParameters{Float64}(
    year_anom = 365.259636 * 86400.0,
    day = 86400.0,
    orbit_semimaj = 1.496e11,
    eccentricity_epoch = 0.0167,
    obliq_epoch = deg2rad(23.44),
    lon_perihelion_epoch = deg2rad(282.95),
    tot_solar_irrad = 1362.0,
    epoch = DateTime(2000, 1, 1, 12, 0, 0),
    mean_anom_epoch = deg2rad(357.5291)
)
```

## Orbital Parameters

Functions for managing time-varying orbital parameters used in paleoclimate studies.

```@docs
Insolation.orbital_params
Insolation.OrbitalDataSplines
```

## Insolation Calculations

Main functions for computing top-of-atmosphere solar radiation.

### Instantaneous Insolation

```@docs
Insolation.insolation
```

### Daily-Averaged Insolation

```@docs
Insolation.daily_insolation
```

### Solar Flux

```@docs
Insolation.solar_flux
```

## Solar Geometry

Functions for calculating solar geometry (distance and solar position in sky). These are typically used internally but can be called directly for specialized applications.

### Instantaneous Geometry

```@docs
Insolation.solar_geometry
```

### Daily-Averaged Geometry

```@docs
Insolation.daily_distance_zenith_angle
```

## Internal Functions

The following functions are typically used internally but are documented for advanced users and developers who need lower-level access to the calculations.

### Temporal and Angular Calculations

```@docs
Insolation.mean_anomaly
Insolation.true_anomaly
Insolation.solar_longitude
Insolation.hour_angle
Insolation.equation_of_time
```

### Distance and Utility Functions

```@docs
Insolation.planet_star_distance
Insolation.years_since_epoch
Insolation.distance_declination_mean_anomaly
```

## Type Aliases

```@docs
Insolation.Parameters.AbstractInsolationParams
```

## Extensions

### CreateParametersExt

The `CreateParametersExt` extension provides integration with [ClimaParams.jl](https://github.com/CliMA/ClimaParams.jl). When both `Insolation.jl` and `ClimaParams.jl` are loaded, this extension enables convenient parameter creation from TOML configuration files, automatically mapping ClimaParams names to Insolation.jl field names.

See the extension source code in `ext/CreateParametersExt.jl` for implementation details.

## Index

```@index
Pages = ["library.md"]
```

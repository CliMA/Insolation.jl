<div align="center">
  <img src="docs/src/assets/logo.svg" alt="Insolation.jl Logo" width="128" height="128">
</div>

# Insolation.jl

|||
|---------------------:|:--------------------------------------------------------------------|
| **Stable Release**   | [![version][version-img]][version-url]                              |
| **Docs Build**       | [![docs build][docs-bld-img]][docs-bld-url]                         |
| **Documentation**    | [![dev][docs-dev-img]][docs-dev-url]                                |
| **Unit Tests**       | [![unit tests][unit-tests-img]][unit-tests-url]                     |
| **Code Coverage**    | [![codecov][codecov-img]][codecov-url]                              |
| **Downloads**        | [![downloads][downloads-img]][downloads-url]                        |

Insolation.jl is a Julia package for calculating solar radiation and solar geometry. It provides efficient and accurate computations of solar zenith angle, azimuth angle, and incoming solar radiation (insolation) for climate science and Earth system modeling applications.

## Features

- **Solar Geometry Calculations**: Compute solar zenith angle, azimuth angle, and declination for any location and time
- **Insolation Computations**: Calculate instantaneous and daily-averaged incoming solar radiation at the top of atmosphere
- **Orbital Parameters**: Support for time-varying orbital parameters using the [Laskar et al. (2004)](https://doi.org/10.1051/0004-6361:20041335) solution for paleoclimate applications
- **Planetary Flexibility**: Works for Earth by default, but applicable to any planetary body with configurable orbital and physical parameters
- **Type-Generic**: Supports different number types (Float32, Float64) for performance optimization
- **Climate Model Integration**: Built for use with [ClimaAtmos.jl](https://github.com/CliMA/ClimaAtmos.jl), [RRTMGP.jl](https://github.com/CliMA/RRTMGP.jl), and other CliMA packages

## Installation

Insolation.jl is a registered Julia package. To install it, open the Julia REPL and run:

```julia
using Pkg
Pkg.add("Insolation")
```

## Quick Start

```julia
using Insolation
using Dates
using ClimaParams

# Calculate insolation for a specific location and time
FT = Float64
lat = FT(40.0)    # degrees North
lon = FT(-105.0)  # degrees East
date = DateTime(2024, 6, 21, 12, 0, 0)  # Summer solstice, noon

# Create Earth parameters with default values
params = InsolationParameters(FT)

# Calculate instantaneous insolation with solar geometry
(; F, S, μ, ζ) = insolation(date, lat, lon, params)
# F: TOA insolation [W m⁻²]
# S: Solar flux [W m⁻²]  
# μ: Cosine of solar zenith angle
# ζ: Solar azimuth angle [radians]

# Calculate daily-averaged insolation
(; F, S, μ) = daily_insolation(date, lat, params)
F_daily, S_daily, μ_daily = F, S, μ
```

For more detailed examples and API documentation, see the [documentation](https://clima.github.io/Insolation.jl/dev/).

[version-img]: https://juliahub.com/docs/General/Insolation/stable/version.svg
[version-url]: https://juliahub.com/ui/Packages/General/Insolation

[docs-bld-img]: https://github.com/CliMA/Insolation.jl/workflows/Documentation/badge.svg
[docs-bld-url]: https://github.com/CliMA/Insolation.jl/actions?query=workflow%3ADocumentation

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://clima.github.io/Insolation.jl/dev/

[unit-tests-img]: https://github.com/CliMA/Insolation.jl/workflows/ci/badge.svg
[unit-tests-url]: https://github.com/CliMA/Insolation.jl/actions?query=workflow%3Aci

[codecov-img]: https://codecov.io/gh/CliMA/Insolation.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/Insolation.jl

[downloads-img]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FInsolation&query=total_requests&label=Downloads
[downloads-url]: https://juliapkgstats.com/pkg/Insolation

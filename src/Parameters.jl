module Parameters

"""
    AbstractInsolationParams

Abstract base type for insolation parameter sets.

All parameter structs used in `Insolation.jl` should inherit from this type.
The main concrete implementation is `InsolationParameters`.

This type hierarchy allows for flexible parameter management and enables
package extensions to provide alternative parameter implementations while
maintaining API compatibility.
"""
abstract type AbstractInsolationParams end
const AIP = AbstractInsolationParams

import Dates: DateTime

"""
    InsolationParameters{FT}

The orbital, solar, and epoch parameters needed for insolation calculations.

This is the main concrete subtype of [`AbstractInsolationParams`](@ref). It can be
constructed directly with the keyword constructor (see the fields below), or, when
`ClimaParams.jl` is loaded, with the convenience constructor `InsolationParameters(FT)`,
which fills in default values for modern Earth.

The values stored here describe the planet's orbit at a fixed reference epoch (typically
J2000); time-varying (Milankovitch) orbital parameters are supplied separately through
[`Insolation.OrbitalDataSplines`](@ref).

# Fields
- `year_anom`: Anomalistic year (perihelion to perihelion) [seconds].
- `day`: Length of a solar day [seconds].
- `eccentricity_epoch`: Eccentricity at epoch [unitless].
- `obliq_epoch`: Obliquity (axial tilt) at epoch [radians].
- `lon_perihelion_epoch`: Longitude of perihelion at epoch [radians].
- `tot_solar_irrad`: Total solar irradiance at the mean orbital distance (semi-major axis) [W m⁻²].
- `epoch`: Reference epoch time [`DateTime`].
- `mean_anom_epoch`: Mean anomaly at epoch [radians].
"""
Base.@kwdef struct InsolationParameters{FT} <: AbstractInsolationParams
    # Orbital periods
    "Anomalistic year (perihelion to perihelion) [seconds]"
    year_anom::FT
    "Length of a solar day [seconds]"
    day::FT

    # Orbital geometry
    "Eccentricity at epoch [unitless]"
    eccentricity_epoch::FT
    "Obliquity at epoch [radians]"
    obliq_epoch::FT
    "Longitude of perihelion (geocentric longitude of the Sun at perihelion relative to vernal equinox) at epoch [radians]"
    lon_perihelion_epoch::FT

    # Solar and Epoch parameters
    "Total solar irradiance at the mean orbital distance (semi-major axis) [W m⁻²]"
    tot_solar_irrad::FT
    "Reference epoch time [DateTime]"
    epoch::DateTime
    "Mean anomaly at epoch [radians]"
    mean_anom_epoch::FT
end

# Make InsolationParameters behave as a scalar in broadcasting operations
# This allows users to write: insolation.(dates, lats, lons, params)
# instead of: insolation.(dates, lats, lons, Ref(params))
Base.broadcastable(x::InsolationParameters) = tuple(x)

# Allows inference of parameter type 
Base.eltype(::InsolationParameters{FT}) where {FT} = FT

# Method wrappers
# This loop creates getter functions for each field, e.g.:
# `year_anom(ps::AIP) = ps.year_anom`
# This allows the rest of the code to use `IP.year_anom(params)`
# instead of `params.year_anom`.
for var in fieldnames(InsolationParameters)
    @eval $var(ps::AIP) = ps.$var
end

end # module

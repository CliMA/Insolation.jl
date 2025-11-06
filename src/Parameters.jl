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

Base.@kwdef struct InsolationParameters{FT} <: AbstractInsolationParams
    # Orbital periods
    "Anomalistic year (perihelion to perihelion) [seconds]"
    year_anom::FT
    "Length of a solar day [seconds]"
    day::FT

    # Orbital geometry
    "Orbit semi-major axis (mean planet-star distance) [m]"
    orbit_semimaj::FT
    "Eccentricity at epoch [unitless]"
    eccentricity_epoch::FT
    "Obliquity at epoch [radians]"
    obliq_epoch::FT
    "Longitude of perihelion (geocentric longitude of the Sun at perihelion relative to vernal equinox) at epoch [radians]"
    lon_perihelion_epoch::FT

    # Solar and Epoch parameters
    "Total Solar Irradiance at 1 au [W m⁻²]"
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

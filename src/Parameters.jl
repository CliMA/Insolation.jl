module Parameters

abstract type AbstractInsolationParams end
const AIP = AbstractInsolationParams

"""
    InsolationParameters

A struct to hold orbital and solar parameters for insolation calculations.
These parameters are typically fixed for a specific epoch (e.g., J2000).
"""
Base.@kwdef struct InsolationParameters{FT, S <: AbstractString} <: AbstractInsolationParams
    # Orbital periods
    "Anomalistic year (perihelion to perihelion) [seconds]"
    year_anom::FT
    "Length of a solar day [seconds]"
    day::FT
    
    # Orbital geometry
    "Orbit semi-major axis (mean Earth-Sun distance) [m]"
    orbit_semimaj::FT
    "Eccentricity at epoch [unitless]"
    eccentricity_epoch::FT
    "Obliquity at epoch [radians]"
    obliq_epoch::FT
    "Heliocentric longitude of perihelion at epoch [radians]"
    lon_perihelion_epoch::FT

    # Solar and Epoch parameters
    "Total Solar Irradiance at 1 au [W m⁻²]"
    tot_solar_irrad::FT
     "Reference epoch time string (e.g., \"2000-01-01T12:00:00.0\")"
    epoch::String
    "Mean anomaly at epoch [radians]"
    mean_anom_epoch::FT
end

# Method wrappers
# This loop creates getter functions for each field, e.g.:
# `year_anom(ps::AIP) = ps.year_anom`
# This allows the rest of the code to use `IP.year_anom(params)`
# instead of `params.year_anom`.for var in fieldnames(InsolationParameters)
    @eval $var(ps::AIP) = ps.$var
end

end # module

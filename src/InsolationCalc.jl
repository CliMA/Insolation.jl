export insolation

"""
    insolation(θ::FT, d::FT, param_set::APS) where {FT <: Real}

returns the insolation given the zenith angle and earth-sun distance
param_set is an AbstractParameterSet from CLIMAParameters.jl
"""
function insolation(θ::FT, d::FT, param_set::APS) where {FT <: Real}
    S0::FT = tot_solar_irrad(param_set)
    AU::FT = astro_unit()

    # set max. zenith angle to π/2, insolation should not be negative
    if θ > π/2
        θ = π/2
    end
    # weighted irradiance (3.12)
    S = S0 * (AU / d)^2
    # TOA insolation (3.15)
    F = S * cos(θ)
    return F
end
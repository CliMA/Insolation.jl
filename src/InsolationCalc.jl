module InsolationCalc

using Dates
using ..ZenithAngleCalc
using ..OrbitalParameters

export insolation

"""
    insolation(θ::FT, d::FT) where {FT <: Real}

returns the insolation given the zenith angle and earth-sun distance
"""
function insolation(θ::FT, d::FT) where {FT <: Real}
    # TOA radiative flux (3.12)
    S = solar_insolation() * (astro_unit() / d)^2
    # set max. zenith angle to π/2, insolation should not be negative
    if θ > π/2
        θ = π/2
    end
    # TOA insolation (3.15)
    F = S * cos(θ)
    return F
end

end
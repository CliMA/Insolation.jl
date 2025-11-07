"""
    CreateParametersExt

Package extension for Insolation.jl that provides integration with ClimaParams.jl.

This extension enables loading insolation parameters from TOML configuration files
using the ClimaParams parameter management system. It provides convenience constructors
for `InsolationParameters` that automatically load parameters from ClimaParams dictionaries,
making it easy to maintain consistent parameter values across CliMA packages.

The extension is automatically loaded when both Insolation.jl and ClimaParams.jl are
present in the environment.
"""
module CreateParametersExt

import Insolation: InsolationParameters

import ClimaParams as CP

"""
    InsolationParameters(::Type{FT}, overrides = NamedTuple()) where {FT <: AbstractFloat}

Create `InsolationParameters` from default ClimaParams TOML configuration.

This convenience constructor loads parameters from the default ClimaParams TOML dictionary
for the specified floating-point type. It automatically maps ClimaParams parameter names
to Insolation.jl field names.

# Arguments
- `FT::Type{<:AbstractFloat}`: The floating-point type for parameter values (e.g., Float64, Float32)
- `overrides::NamedTuple`: Optional named tuple to override specific parameter values

# Returns
- `InsolationParameters{FT}`: Parameter struct containing all orbital and solar parameters

# Example
```julia
using Insolation
using ClimaParams

# Create parameters with default values
params = InsolationParameters(Float64)

# Create parameters with custom solar irradiance
params = InsolationParameters(Float64, (; tot_solar_irrad = 1365.0))
```
"""
InsolationParameters(::Type{FT}, overrides = NamedTuple()) where {FT<:AbstractFloat} =
    InsolationParameters(CP.create_toml_dict(FT), overrides)

"""
    InsolationParameters(toml_dict::CP.ParamDict{FT}, overrides = NamedTuple()) where {FT}

Create `InsolationParameters` from a ClimaParams parameter dictionary.

Constructs an `InsolationParameters` object by extracting and mapping parameters from
a ClimaParams dictionary. This method handles the translation between ClimaParams
naming conventions and Insolation.jl field names.

# Arguments
- `toml_dict::CP.ParamDict{FT}`: ClimaParams parameter dictionary
- `overrides::NamedTuple`: Optional named tuple to override specific parameter values

# Returns
- `InsolationParameters{FT}`: Parameter struct containing all orbital and solar parameters

# Parameter Mapping
The following ClimaParams names are mapped to Insolation.jl fields:
- `epoch_time` → `epoch`: Reference epoch (DateTime)
- `day` → `day`: Length of solar day [s]
- `anomalistic_year_length` → `year_anom`: Anomalistic year length [s]
- `length_orbit_semi_major` → `orbit_semimaj`: Orbital semi-major axis [m]
- `orbit_eccentricity_at_epoch` → `eccentricity_epoch`: Eccentricity at epoch
- `total_solar_irradiance` → `tot_solar_irrad`: Solar irradiance at 1 AU [W m⁻²]
- `orbit_obliquity_at_epoch` → `obliq_epoch`: Obliquity at epoch [rad]
- `mean_anomaly_at_epoch` → `mean_anom_epoch`: Mean anomaly at epoch [rad]
- `longitude_perihelion_at_epoch` → `lon_perihelion_epoch`: Longitude of perihelion at epoch [rad]

# Example
```julia
using Insolation
using ClimaParams as CP

# Load parameters from custom TOML file
toml_dict = CP.create_toml_dict(Float64; param_file = "my_params.toml")
params = InsolationParameters(toml_dict)

# Override specific parameters
params = InsolationParameters(toml_dict, (; tot_solar_irrad = 1362.0))
```
"""
function InsolationParameters(
    toml_dict::CP.ParamDict{FT},
    overrides = NamedTuple(),
) where {FT}
    name_map = (;
        :epoch_time => :epoch,
        :day => :day,
        :anomalistic_year_length => :year_anom,
        :length_orbit_semi_major => :orbit_semimaj,
        :orbit_eccentricity_at_epoch => :eccentricity_epoch,
        :total_solar_irradiance => :tot_solar_irrad,
        :orbit_obliquity_at_epoch => :obliq_epoch,
        :mean_anomaly_at_epoch => :mean_anom_epoch,
        :longitude_perihelion_at_epoch => :lon_perihelion_epoch,
    )

    parameters = CP.get_parameter_values(toml_dict, name_map, "Insolation")
    parameters = merge(parameters, overrides)
    return InsolationParameters{FT}(; parameters...)
end

end

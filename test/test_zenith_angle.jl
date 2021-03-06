using Dates
using Insolation 

using CLIMAParameters
using CLIMAParameters.Planet
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

rtol = 1e-2

## Test Zenith
date = Dates.DateTime(2020, 1, 1, 6, 0, 0)
lon, lat = [0.0, 0.0]
sza, azi, d = instantaneous_zenith_angle(date, lon, lat, param_set)
@test sza ≈ π/2 rtol=rtol

date = Dates.DateTime(2020, 1, 1, 18, 0, 0)
sza, azi, d = instantaneous_zenith_angle(date, lon, lat, param_set)
@test sza ≈ π/2 rtol=rtol

## Test Azimuth
date = Dates.DateTime(2020, 1, 1, 0, 0, 0)
lon, lat = [0.0, 40.0]
sza, azi, d = instantaneous_zenith_angle(date, lon, lat, param_set)
@test azi ≈ π/2 rtol=rtol

date = Dates.DateTime(2020, 1, 1, 12, 0, 0)
sza, azi, d = instantaneous_zenith_angle(date, lon, lat, param_set)
@test azi ≈ 3π/2 rtol=rtol

## Test Distance
date = Dates.DateTime(2000, 3, 22, 0, 0, 0)
sza, azi, d = instantaneous_zenith_angle(date, lon, lat, param_set)
@test d ≈ orbit_semimaj(param_set) rtol=rtol
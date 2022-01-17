push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Dates
using Insolation 

using CLIMAParameters
using CLIMAParameters.Planet
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

rtol = 1e-2

# sunrise at equator
date = Dates.DateTime(2020, 2, 20, 6, 11, 0)
lon, lat = [0.0, 0.0]
sza, azi, d = instantaneous_zenith_angle(date, lon, lat, param_set)
@test sza ≈ π/2 rtol=rtol

# solar noon at equator
date = Dates.DateTime(2020, 2, 20, 12, 13, 0)
sza, azi, d = instantaneous_zenith_angle(date, lon, lat, param_set)
@test azi ≈ 3π/2 rtol=rtol

# sunset at equator
date = Dates.DateTime(2020, 2, 20, 18, 17, 0)
sza, azi, d = instantaneous_zenith_angle(date, lon, lat, param_set)
@test sza ≈ π/2 rtol=rtol

# solar noon at equator, eot correction = false
date = Dates.DateTime(2020, 2, 20, 12, 0, 0)
sza, azi, d = instantaneous_zenith_angle(date, lon, lat, param_set; eot_correction=false)
println(sza, azi)
@test azi ≈ 3π/2 rtol=rtol

# sunset at equator, eot correction = false
date = Dates.DateTime(2020, 2, 20, 18, 0, 0)
sza, azi, d = instantaneous_zenith_angle(date, lon, lat, param_set; eot_correction=false)
println(sza, azi)
@test sza ≈ π/2 rtol=rtol

## Test Polar Night
# polar night NH 1
date = Dates.DateTime(2020, 12, 20, 11, 0, 0)
lon, lat = [0.0, 80.0]
sza, azi, d = instantaneous_zenith_angle(date, lon, lat, param_set)
@test sza > π/2

# polar night NH 2
date = Dates.DateTime(2020, 12, 20, 23, 0, 0)
lon, lat = [0.0, 80.0]
sza, azi, d = instantaneous_zenith_angle(date, lon, lat, param_set)
@test sza > π/2

# polar night SH 1
date = Dates.DateTime(2020, 6, 20, 11, 0, 0)
lon, lat = [0.0, -80.0]
sza, azi, d = instantaneous_zenith_angle(date, lon, lat, param_set)
@test sza > π/2

# polar night SH 2
date = Dates.DateTime(2020, 6, 20, 23, 0, 0)
lon, lat = [0.0, -80.0]
sza, azi, d = instantaneous_zenith_angle(date, lon, lat, param_set)
@test sza > π/2

## Test Distance
date = Dates.DateTime(2000, 3, 22, 0, 0, 0)
sza, azi, d = instantaneous_zenith_angle(date, lon, lat, param_set)
@test d ≈ orbit_semimaj(param_set) rtol=rtol
date = Dates.DateTime(2020, 2, 20, 11, 11, 0)
lon, lat = [FT(80.0), FT(20.0)]

od = Insolation.OrbitalData()

sza, azi, d = instantaneous_zenith_angle(
    date,
    date0,
    lon,
    lat,
    param_set,
    eot_correction = false,
)


F = insolation(sza, d, param_set)
@test typeof(sza) == FT
@test typeof(azi) == FT
@test typeof(d) == FT
@test typeof(F) == FT

sza, azi, d = instantaneous_zenith_angle(
    date,
    date0,
    od,
    lon,
    lat,
    param_set,
    eot_correction = false,
)

F = insolation(sza, d, param_set)
@test typeof(sza) == FT
@test typeof(azi) == FT
@test typeof(d) == FT
@test typeof(F) == FT

sza, azi, d = instantaneous_zenith_angle(
    date,
    date0,
    lon,
    lat,
    param_set,
    eot_correction = true,
)

F = insolation(sza, d, param_set)
@test typeof(sza) == FT
@test typeof(azi) == FT
@test typeof(d) == FT
@test typeof(F) == FT

sza, azi, d = instantaneous_zenith_angle(
    date,
    date0,
    od,
    lon,
    lat,
    param_set,
    eot_correction = true,
)


F = insolation(sza, d, param_set)
@test typeof(sza) == FT
@test typeof(azi) == FT
@test typeof(d) == FT
@test typeof(F) == FT

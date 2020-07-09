using Dates
using Insolation.SolarZenithAngle

rtol = 1e-2

PST = -7.0
California = [-105.0, 40.0]

tz = PST
lon, lat = California

date = Dates.DateTime(2020,3,21,12,0,0)
sza = instantaneous_zenith_angle(date, tz, lon, lat)
@test rad2deg(sza) ≈ 20.0 rtol=rtol

γ = 23.44
ϖ = 282.95
e = 0.017
sza = instantaneous_zenith_angle(date, tz, lon, lat, γ, ϖ, e)
@test rad2deg(sza) ≈ 20.0 rtol=rtol

sza = daily_zenith_angle(date, lat)
@test rad2deg(sza) ≈ 20.0 rtol=rtol

days_since_equinox = 85
sza = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
@test rad2deg(sza) ≈ 20.0 rtol=rtol

ϖ = 282.95 + 180.0
sza = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
@test rad2deg(sza) ≈ 20.0 rtol=rtol

γ = 22.0
ϖ = 282.95
sza = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
@test rad2deg(sza) ≈ 20.0 rtol=rtol

γ = 18.0
sza = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
@test rad2deg(sza) ≈ 20.0 rtol=rtol

γ = 60.0
sza = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
@test rad2deg(sza) ≈ 20.0 rtol=rtol

γ = 97.86
sza = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
@test rad2deg(sza) ≈ 20.0 rtol=rtol
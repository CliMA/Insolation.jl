using Dates
using Insolation.SolarZenithAngle

rtol = 1e-2

PST = -7.0
California = [-105.0, 40.0]

tz = PST
lon, lat = California

date = Dates.DateTime(2020,3,21,12,0,0)

γ = 23.44
ϖ = 282.95
e = 0.017

days_since_equinox = 85

sza, d = instantaneous_zenith_angle(date, tz, lon, lat)
@test rad2deg(sza) ≈ 39.573 rtol=rtol

sza, d = instantaneous_zenith_angle(date, tz, lon, lat, γ, ϖ, e)
@test rad2deg(sza) ≈ 76.718 rtol=rtol

sza, d = daily_zenith_angle(date, lat)
@test rad2deg(sza) ≈ 75.737 rtol=rtol

sza, d = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
@test rad2deg(sza) ≈ 68.559 rtol=rtol

ϖ = 282.95 + 180.0
sza = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
@test rad2deg(sza) ≈ 68.982 rtol=rtol

γ = 22.0
ϖ = 282.95
sza, d = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
println(rad2deg(sza))
@test rad2deg(sza) ≈ 68.982 rtol=rtol

γ = 18.0
sza, d = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
println(rad2deg(sza))
@test rad2deg(sza) ≈ 70.177 rtol=rtol

γ = 60.0
sza, d = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
println(rad2deg(sza))
@test rad2deg(sza) ≈ 56.495 rtol=rtol

γ = 97.86
sza, d = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
println(rad2deg(sza))
@test rad2deg(sza) ≈ 50.846 rtol=rtol
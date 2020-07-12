using Dates
using Insolation.SolarZenithAngle
include("reference_sza_codes.jl")

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

println("instantaneous tests")
sza, d = instantaneous_zenith_angle(date, tz, lon, lat)
sza_ref = inst_sza(date, tz, lon, lat)
println(rad2deg(sza), "\t", sza_ref)
@test rad2deg(sza) ≈ 39.574 rtol=rtol

date = Dates.DateTime(2020,3,21,17,0,0)
sza, d = instantaneous_zenith_angle(date, tz, lon, lat)
sza_ref = inst_sza(date, tz, lon, lat)
println(rad2deg(sza), "\t", sza_ref)
@test rad2deg(sza) ≈ 77.037 rtol=rtol

date = Dates.DateTime(2020,3,21,12,0,0)
lat = 85.0
sza, d = instantaneous_zenith_angle(date, tz, lon, lat)
sza_ref = inst_sza(date, tz, lon, lat)
println(rad2deg(sza), "\t", sza_ref)
@test rad2deg(sza) ≈ 84.551 rtol=rtol

date = Dates.DateTime(2020,3,21,4,0,0)
lon, lat = California
sza, d = instantaneous_zenith_angle(date, tz, lon, lat)
# sza_ref = inst_sza(date, tz, lon, lat)
# println(rad2deg(sza), "\t", sza_ref)
@test rad2deg(sza) ≈ 90.000 rtol=rtol

date = Dates.DateTime(2020,8,21,11,0,0)
lat = -85.0
sza, d = instantaneous_zenith_angle(date, tz, lon, lat)
# sza_ref = inst_sza(date, tz, lon, lat)
# println(rad2deg(sza), "\t", sza_ref)
@test rad2deg(sza) ≈ 90.000 rtol=rtol

lat = 40.0
lon = 100.0
sza, d = instantaneous_zenith_angle(date, tz, lon, lat)
sza_ref = inst_sza(date, tz, lon, lat)
println(rad2deg(sza), "\t", sza_ref)
#@test rad2deg(sza) ≈ 90.000 rtol=rtol

date = Dates.DateTime(2020,3,21,12,0,0)
lon, lat = California
sza, d = instantaneous_zenith_angle(date, tz, lon, lat, γ, ϖ, e)
@test rad2deg(sza) ≈ 76.718 rtol=rtol

sza, d = daily_zenith_angle(date, lat)
@test rad2deg(sza) ≈ 75.737 rtol=rtol

lon, lat = California
sza, d = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
# sza_ref = daily_sza(days_since_equinox, lat, γ, ϖ, e)
# println(rad2deg(sza), "\t", sza_ref)
@test rad2deg(sza) ≈ 68.559 rtol=rtol

lat = -85.0
sza, d = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
# sza_ref = daily_sza(days_since_equinox, lat, γ, ϖ, e)
# println(rad2deg(sza), "\t", sza_ref)
@test rad2deg(sza) ≈ 90.000 rtol=rtol

lon, lat = California
ϖ = 282.95 + 180.0
sza, d = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
# sza_ref = daily_sza(days_since_equinox, lat, γ, ϖ, e)
# println(rad2deg(sza), "\t", sza_ref)
@test rad2deg(sza) ≈ 68.524 rtol=rtol

γ = 22.0
ϖ = 282.95
sza, d = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
# sza_ref = daily_sza(days_since_equinox, lat, γ, ϖ, e)
# println(rad2deg(sza), "\t", sza_ref)
@test rad2deg(sza) ≈ 68.982 rtol=rtol

γ = 18.0
sza, d = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
# sza_ref = daily_sza(days_since_equinox, lat, γ, ϖ, e)
# println(rad2deg(sza), "\t", sza_ref)
@test rad2deg(sza) ≈ 70.177 rtol=rtol

γ = 60.0
sza, d = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
# sza_ref = daily_sza(days_since_equinox, lat, γ, ϖ, e)
# println(rad2deg(sza), "\t", sza_ref)
@test rad2deg(sza) ≈ 56.495 rtol=rtol

γ = 97.86
sza, d = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
# sza_ref = daily_sza(days_since_equinox, lat, γ, ϖ, e)
# println(rad2deg(sza), "\t", sza_ref)
@test rad2deg(sza) ≈ 50.846 rtol=rtol

# julia --project --color=yes --check-bounds=yes -e 'using Pkg; Pkg.instantiate(); Pkg.build(); Pkg.test(coverage=true);'
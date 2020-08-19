using Dates
using Insolation.ZenithAngleCalc
include("reference_sza_codes.jl")

rtol = 1e-2

PST = -7.0
California = [-105.0, 40.0]

tz = PST
lon, lat = California

γ = 23.44
ϖ = 282.95
e = 0.017

days_since_equinox = 85

doprint = false

date = Dates.DateTime(2020,3,21,12,0,0)
sza, d = instantaneous_zenith_angle(date, tz, lon, lat)
sza_ref = inst_sza(date, tz, lon, lat)
@test rad2deg(sza) ≈ sza_ref rtol=rtol
if doprint
    println(rad2deg(sza), "\t", sza_ref)
end

date = Dates.DateTime(2020,3,21,17,0,0)
sza, d = instantaneous_zenith_angle(date, tz, lon, lat)
sza_ref = inst_sza(date, tz, lon, lat)
@test rad2deg(sza) ≈ sza_ref rtol=rtol
if doprint
    println(rad2deg(sza), "\t", sza_ref)
end

date = Dates.DateTime(2020,3,21,12,0,0)
lat = 85.0
sza, d = instantaneous_zenith_angle(date, tz, lon, lat)
sza_ref = inst_sza(date, tz, lon, lat)
@test rad2deg(sza) ≈ sza_ref rtol=rtol
if doprint
    println(rad2deg(sza), "\t", sza_ref)
end

date = Dates.DateTime(2020,3,21,4,0,0)
lon, lat = California
sza, d = instantaneous_zenith_angle(date, tz, lon, lat)
sza_ref = inst_sza(date, tz, lon, lat)
@test rad2deg(sza) ≈ sza_ref rtol=rtol
if doprint
    println(rad2deg(sza), "\t", sza_ref)
end

date = Dates.DateTime(2020,8,21,11,0,0)
lat = -85.0
sza, d = instantaneous_zenith_angle(date, tz, lon, lat)
sza_ref = inst_sza(date, tz, lon, lat)
@test rad2deg(sza) ≈ sza_ref rtol=rtol
if doprint
    println(rad2deg(sza), "\t", sza_ref)
end

date = Dates.DateTime(2020,8,21,11,0,0)
lon, lat = [105.0, 40.0]
tz = 7.0
sza, d = instantaneous_zenith_angle(date, tz, lon, lat)
sza_ref = inst_sza(date, tz, lon, lat)
@test rad2deg(sza) ≈ sza_ref rtol=rtol
if doprint
    println(rad2deg(sza), "\t", sza_ref)
end

date = Dates.DateTime(2020,8,21,11,0,0)
lon, lat = [-105.0, 40.0]
tz = -7.0
sza, d = instantaneous_zenith_angle(date, tz, lon, lat)
sza_ref = inst_sza(date, tz, lon, lat)
@test rad2deg(sza) ≈ sza_ref rtol=rtol
if doprint
    println(rad2deg(sza), "\t", sza_ref)
end

date = Dates.DateTime(2020,3,21,12,0,0)
lon, lat = California
tz = PST
sza, d = instantaneous_zenith_angle(date, tz, lon, lat, γ, ϖ, e)
sza_ref = inst_sza(date, tz, lon, lat)
@test rad2deg(sza) ≈ sza_ref rtol=rtol
if doprint
    println(rad2deg(sza), "\t", sza_ref)
end

sza, d = daily_zenith_angle(date, lat)
sza_ref = daily_sza(1, lat, γ, ϖ, e)
@test rad2deg(sza) ≈ sza_ref rtol=rtol
if doprint
    println(rad2deg(sza), "\t", sza_ref)
end

lon, lat = California
sza, d = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
sza_ref = daily_sza(days_since_equinox, lat, γ, ϖ, e)
@test rad2deg(sza) ≈ sza_ref rtol=rtol
if doprint
    println(rad2deg(sza), "\t", sza_ref)
end

lat = -85.0
sza, d = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
sza_ref = daily_sza(days_since_equinox, lat, γ, ϖ, e)
@test rad2deg(sza) ≈ sza_ref rtol=rtol
if doprint
    println(rad2deg(sza), "\t", sza_ref)
end

lon, lat = California
ϖ = 282.95 + 180.0
sza, d = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
sza_ref = daily_sza(days_since_equinox, lat, γ, ϖ, e)
@test rad2deg(sza) ≈ sza_ref rtol=rtol
if doprint
    println(rad2deg(sza), "\t", sza_ref)
end

γ = 22.0
ϖ = 282.95
sza, d = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
sza_ref = daily_sza(days_since_equinox, lat, γ, ϖ, e)
@test rad2deg(sza) ≈ sza_ref rtol=rtol
if doprint
    println(rad2deg(sza), "\t", sza_ref)
end

γ = 18.0
sza, d = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
sza_ref = daily_sza(days_since_equinox, lat, γ, ϖ, e)
@test rad2deg(sza) ≈ sza_ref rtol=rtol
if doprint
    println(rad2deg(sza), "\t", sza_ref)
end

γ = 97.86
sza, d = daily_zenith_angle(days_since_equinox, γ, ϖ, e, lat)
sza_ref = daily_sza(days_since_equinox, lat, γ, ϖ, e)
@test rad2deg(sza) ≈ sza_ref rtol=rtol
if doprint
    println(rad2deg(sza), "\t", sza_ref)
end

# julia --project --color=yes --check-bounds=yes -e 'using Pkg; Pkg.instantiate(); Pkg.build(); Pkg.test(coverage=true);'
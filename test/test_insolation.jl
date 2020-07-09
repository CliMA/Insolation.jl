using Dates
using Insolation.SolarInsolation

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

S = instantaneous_insolation(date, tz, lon, lat)
@test S ≈ 1057.824 rtol=rtol

S = daily_insolation(date, lat)
@test S ≈ 338.117 rtol=rtol

θ = 58.0
d = 1.50e11 
S = insolation(θ, d)
@test S ≈ 161.454 rtol=rtol

S = daily_insolation(days_since_equinox, γ, ϖ, e, lat)
@test S ≈ 482.404 rtol=rtol

ϖ = 282.95 + 180.0
S = daily_insolation(days_since_equinox, γ, ϖ, e, lat)
@test S ≈ 515.197 rtol=rtol

γ = 22.0
ϖ = 282.95
S = daily_insolation(days_since_equinox, γ, ϖ, e, lat)
@test S ≈ 473.324 rtol=rtol

γ = 18.0
S = daily_insolation(days_since_equinox, γ, ϖ, e, lat)
@test S ≈ 447.527 rtol=rtol

γ = 60.0
S = daily_insolation(days_since_equinox, γ, ϖ, e, lat)
@test S ≈ 728.479 rtol=rtol

γ = 97.86
S = daily_insolation(days_since_equinox, γ, ϖ, e, lat)
@test S ≈ 833.273 rtol=rtol
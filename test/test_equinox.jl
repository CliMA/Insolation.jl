push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Statistics
using Roots
using Dates
using Insolation 
using Plots

using CLIMAParameters
using CLIMAParameters.Planet
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()
#CLIMAParameters.Planet.lon_perihelion(::EarthParameterSet) = deg2rad(282.937348 + 180)

atol = 1e-6
rtol = 1e-2

# Difference in NH and SH zenith angles at time x in given year
function zdiff(x, year)
    date = xtodate(x,year)
    theta_s, dist = daily_zenith_angle(date, -45., param_set)
    theta_n, dist = daily_zenith_angle(date, 45., param_set)
    return theta_n - theta_s
end

# x is date relative to March 1, with 1.00 representing March 1 00:00
function xtodate(x, year)
    basedate = Dates.DateTime(year, 3, 1)
    deltat = Dates.Second(round((x-1)*86400))
    date = basedate + deltat
    return date
end

day = zeros(length(1900:2100))
for (i,year) in enumerate(1900:2100)
    f = (x -> zdiff(x, year))
    # println(f(1), f(15), f(30))
    day[i] = find_zeros(f,1.,30)[1]
end

# plot(1900:2100, day)
# xlabel!("Year")
# ylabel!("Day in March")
# title!("Date of vernal equinox")

# test mean is about March 21
@test mean(day) ≈ 21 atol=0.5

test decreasing
@test mean(day[:100]) > mean(day[100:end])

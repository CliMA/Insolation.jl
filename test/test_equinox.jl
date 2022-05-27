push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Statistics
using Roots
using Dates
using Insolation 
using Plots

using CLIMAParameters
using CLIMAParameters.Planet: year_anom
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()
CLIMAParameters.Planet.year_anom(::EarthParameterSet) = 365.24219 * CLIMAParameters.Planet.day(param_set)
# 365.256363004 = sidereal, 365.24219 tropical, 365.259636 anomalistic

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
    return basedate + deltat
end

days = zeros(length(1900:2100))
for (i,year) in enumerate(1900:2100)
    f = (x -> zdiff(x, year))
    days[i] = find_zeros(f,1.,30)[1]
end

# test mean is about March 21
@test mean(days) ≈ 21 atol=0.5

# test decreasing
@test mean(days[:100]) > mean(days[100:end])

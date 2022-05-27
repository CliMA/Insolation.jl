push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Statistics
using Optim
using Dates
using Insolation 
using Plots

using CLIMAParameters
using CLIMAParameters.Planet
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()
CLIMAParameters.Planet.year_anom(::EarthParameterSet) = 365.24219 * CLIMAParameters.Planet.day(param_set)
# 365.256363004 = sidereal, 365.24219 tropical, 365.259636 anomalistic

atol = 1e-6
rtol = 1e-2

# x is date relative to Jan 1, with 1.00 representing Jan 1 00:00
function xtojandate(x, year)
    basedate = Dates.DateTime(year, 1, 1)
    deltat = Dates.Second(round((x-1)*86400))
    date = basedate + deltat
    return date
end

# Earth-Sun distance
function edist(x, year)
    date = xtojandate(x,year)
    _, dist = daily_zenith_angle(date, 0., param_set)
    return dist/astro_unit()
end

years = 1900:2100
days = zeros(length(years))
for (i,year) in enumerate(years)
    f = (x -> edist(x, year))
    res = optimize(f,1.,30)
    days[i] = Optim.minimizer(res)[1]
end

# test mean is about Jan 3.5
@test mean(days) ≈ 3.5 atol=0.5

# test increasing
@test mean(days[:100]) < mean(days[100:end])

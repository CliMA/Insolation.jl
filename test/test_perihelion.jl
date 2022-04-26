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
#CLIMAParameters.Planet.lon_perihelion(::EarthParameterSet) = deg2rad(282.937348 + 180)

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
    theta, dist = daily_zenith_angle(date, 0., param_set)
    return dist/astro_unit()
end

years = 1900:2100
day = zeros(length(years))
for (i,year) in enumerate(years)
    f = (x -> edist(x, year))
    res = optimize(f,1.,30)
    #println(f(1), f(20), f(360))
    day[i] = Optim.minimizer(res)[1]
end

# #print(day)
# plot((years), day)
# xlabel!("Year")
# ylabel!("Day in Jan")
# title!("Date of perihelion")

# test mean is about Jan 3.5
@test mean(day) ≈ 3.5 atol=0.5

# test increasing
@test mean(day[:100]) < mean(day[100:end])

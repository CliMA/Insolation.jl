push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Optim
using Dates
using Insolation 
using Plots

using CLIMAParameters
using CLIMAParameters.Planet
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

atol = 1e-6
rtol = 1e-2

function xtojandate(x, year)
    # x is date relative to Jan 1, with 1.00 representing Jan 1 00:00
    basedate = Dates.DateTime(year, 1, 1)
    deltat = Dates.Second(round((x-1)*86400))
    date = basedate + deltat
    return date
end

function edist(x, year)
    date = xtojandate(x,year)
    theta, dist = daily_zenith_angle(date, 0., param_set)
    return dist/astro_unit()
end

years = 1900:2100
day = zeros(length(years))
for (i,year) in enumerate(years)
    # Temporary function to allow passing the year argument
    f = (x -> edist(x, year))
    res = optimize(f,1.,30)
    day[i] = Optim.minimizer(res)[1]
    # @printf("%d  %6.3f  %s\n", year, d[1], Dates.format(xtojandate(d[1], year), "yyyy-mm-dd HH:MM"))
end

plot((years), day)

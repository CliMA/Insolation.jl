push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Roots
using Dates
using Insolation 
using Plots

using CLIMAParameters
using CLIMAParameters.Planet
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

atol = 1e-6
rtol = 1e-2

function zdiff(x, year)
    # Difference in NH and SH zenith angles at time x in given year
    date = xtodate(x,year)
    theta_s, dist = daily_zenith_angle(date, -45., param_set)
    theta_n, dist = daily_zenith_angle(date, 45., param_set)
    return theta_n - theta_s
end

function xtodate(x, year)
    # x is date relative to March 1, with 1.00 representing March 1 00:00
    basedate = Dates.DateTime(year, 3, 1)
    deltat = Dates.Second(round((x-1)*86400))
    date = basedate + deltat
    return date
end

day = zeros(length(1900:2100))
for (i,year) in enumerate(1900:2100)
    # Temporary function to allow passing the year argument
    f = (x -> zdiff(x, year))
    day[i] = find_zeros(f,1.,30)[1]
    #@printf("%d  %6.3f  %s\n", year, d[1], Dates.format(xtodate(d[1], year), "yyyy-mm-dd HH:MM"))
end

plot(1900:2100, day)

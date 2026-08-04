using Plots
using Dates

using Insolation
import Insolation.Parameters as IP
import ClimaParams as CP

FT = Float64
param_set = IP.InsolationParameters(FT)

# Difference in NH and SH zenith angles at time x in given year
function zdiff(x, year, od)
    date = xtomarchdate(x, year)

    # Get orbital parameters for this time
    Δt_years = Insolation.years_since_epoch(param_set, date)
    orb_params = Insolation.orbital_params(od, Δt_years)

    # Calculate zenith angles for Southern and Northern mid-latitudes
    result_s = Insolation.daily_distance_zenith_angle(date, -45.0, orb_params, param_set)
    result_n = Insolation.daily_distance_zenith_angle(date, 45.0, orb_params, param_set)

    return result_n.daily_θ - result_s.daily_θ
end

# x is date relative to March 1, with 1.00 representing March 1 00:00
function xtomarchdate(x, year)
    basedate = Dates.DateTime(year, 3, 1)
    deltat = Dates.Second(round((x - 1) * IP.day(param_set)))
    return basedate + deltat
end

# Planet-star distance
function edist(x, year, od)
    date = xtojandate(x, year)

    # Get orbital parameters for this time
    Δt_years = Insolation.years_since_epoch(param_set, date)
    orb_params = Insolation.orbital_params(od, Δt_years)

    # Calculate distance
    result = Insolation.daily_distance_zenith_angle(date, 0.0, orb_params, param_set)

    return result.d / IP.orbit_semimaj(param_set)
end

# x is date relative to Jan 1, with 1.00 representing Jan 1 00:00
function xtojandate(x, year)
    basedate = Dates.DateTime(year, 1, 1)
    deltat = Dates.Second(round((x - 1) * IP.day(param_set)))
    date = basedate + deltat
    return date
end

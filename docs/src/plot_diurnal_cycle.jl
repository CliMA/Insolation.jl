using Insolation
using Plots
using Dates
using Statistics
using Formatting

using CLIMAParameters
using CLIMAParameters.Planet
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using CLIMAParameters: AbstractParameterSet
const APS = AbstractParameterSet

function diurnal_cycle(lat, lon, date, timezone, filename)
    nhours = 1000
    hours = collect(range(0, stop = 24, length = nhours))
    insol = zeros(nhours)
    sza = zeros(nhours)
    for (i,hr) in enumerate(hours)
        h = Int(round(hr + timezone))
        m = Int(round((hr + timezone - h) * 60))
        datetime = date + Dates.Hour(h) + Dates.Minute(m)
        S, mu = solar_flux_and_cos_sza(datetime, lon, lat, param_set)
        insol[i] = S*mu
        sza[i] = rad2deg(acos(mu))
    end
    plot(hours, insol, color=:blue, lw=2, label="", 
            ylabel="Insolation [W/m2]", yguidefontcolor=:blue, dpi=200, 
            left_margein=10Plots.mm, right_margin=15Plots.mm)
    plot!(twinx(), hours, sza, color=:red, lw=2, label="", 
            ylabel="SZA [deg]", yguidefontcolor=:red)
    xlabel!("Local Time")
    savefig(filename)
end

# Milankovitch Cycles

## Variations in orbital parameters
```@example
using Plots #hide
using Dates #hide
using CLIMAParameters #hide
using CLIMAParameters.Planet #hide
struct EarthParameterSet <: AbstractEarthParameterSet end #hide
const param_set = EarthParameterSet() #hide
CLIMAParameters.Planet.year_anom(::EarthParameterSet) = 365.24219 * CLIMAParameters.Planet.day(param_set)

Ya = year_anom(param_set)
ϖ0 = lon_perihelion_epoch(param_set)
γ0 = obliq_epoch(param_set)
e0 = eccentricity_epoch(param_set)

dt = collect(-300e3:300e3) .* Ya;
ϖ = mod.(ϖ0 .+ 2π*dt/(26e3*Ya) + 2π*dt/(112e3*Ya), 2π);
γ = γ0 .+ deg2rad(1.2)*sin.(2π*dt/(41e3*Ya));
e = e0 .+ 0.02*sin.(2π*dt/(405e3*Ya)) .+ 0.01*sin.(2π*dt/(110e3*Ya)) .+ 0.01;

p1 = plot(dt ./ (1e3*Ya), sin.(ϖ), legend=false);
ylabel!("sin(ϖ)");
p2 = plot(dt ./ (1e3*Ya), γ, legend=false);
ylabel!("γ");
p3 = plot(dt ./ (1e3*Ya), e, legend=false);
ylabel!("e");
xlabel!("time (kY)")
plot(p1, p2, p3, layout = grid(3,1), size=(600,400), dpi=150);
savefig("orbital_params.png")
```
![](orbital_params.png)

## Slow variations in date of vernal equinox and perihelion
```@example
using Insolation #hide
using Plots #hide
using Dates #hide
using Roots #hide
using Optim #hide
using CLIMAParameters #hide
using CLIMAParameters.Planet #hide
struct EarthParameterSet <: AbstractEarthParameterSet end #hide
const param_set = EarthParameterSet() #hide
CLIMAParameters.Planet.year_anom(::EarthParameterSet) = 365.24219 * CLIMAParameters.Planet.day(param_set)

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

plot(1900:2100, days)
xlabel!("Year")
ylabel!("Day in March")
title!("Date of vernal equinox")
savefig("equinox_dates.png")
```
![](equinox_dates.png)

```@example
using Insolation #hide
using Plots #hide
using Dates #hide
using Roots #hide
using Optim #hide
using CLIMAParameters #hide
using CLIMAParameters.Planet #hide
struct EarthParameterSet <: AbstractEarthParameterSet end #hide
const param_set = EarthParameterSet() #hide
CLIMAParameters.Planet.year_anom(::EarthParameterSet) = 365.24219 * CLIMAParameters.Planet.day(param_set)

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

plot((years), days)
xlabel!("Year")
ylabel!("Day in Jan")
title!("Date of perihelion")
savefig("perihelion_dates.png")
```
![](perihelion_dates.png)
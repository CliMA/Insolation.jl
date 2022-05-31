# Milankovitch Cycles

## Variations in orbital parameters
```@example
using Insolation #hide
using Plots #hide
using CLIMAParameters #hide
using CLIMAParameters.Planet #hide
struct EarthParameterSet <: AbstractEarthParameterSet end #hide
const param_set = EarthParameterSet() #hide

dt = collect(-300e3:100:300e3);
y = transpose(reduce(hcat, orbital_params.(dt)))
ϖ, γ, e = y[:,1], y[:,2], y[:,3]

p1 = plot(dt ./ (1e3), sin.(ϖ), legend=false);
ylabel!("sin(ϖ)");
p2 = plot(dt ./ (1e3), γ, legend=false);
ylabel!("γ");
p3 = plot(dt ./ (1e3), e, legend=false);
ylabel!("e");
xlabel!("time (kY)")
plot(p1, p2, p3, layout = grid(3,1), size=(600,400), dpi=150);
savefig("orbital_params.png")
```
![](orbital_params.png)

## Variations in date of vernal equinox and perihelion on centennial timescales
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

mm = false

# Difference in NH and SH zenith angles at time x in given year
function zdiff(x, year)
    date = xtodate(x,year)
    theta_s, dist = daily_zenith_angle(date, -45., param_set, milankovitch=mm)
    theta_n, dist = daily_zenith_angle(date, 45., param_set, milankovitch=mm)
    return theta_n - theta_s
end

# x is date relative to March 1, with 1.00 representing March 1 00:00
function xtodate(x, year)
    basedate = Dates.DateTime(year, 3, 1)
    deltat = Dates.Second(round((x-1)*86400))
    return basedate + deltat
end

# Earth-Sun distance
function edist(x, year)
    date = xtojandate(x,year)
    _, dist = daily_zenith_angle(date, 0., param_set, milankovitch=mm)
    return dist/astro_unit()
end

# x is date relative to Jan 1, with 1.00 representing Jan 1 00:00
function xtojandate(x, year)
    basedate = Dates.DateTime(year, 1, 1)
    deltat = Dates.Second(round((x-1)*86400))
    date = basedate + deltat
    return date
end

years = 1800:2200
days_eq = zeros(length(years))
days_per = zeros(length(years))
for (i,year) in enumerate(years)
    f = (x -> zdiff(x, year))
    days_eq[i] = find_zeros(f,-30,60)[1]

    f = (x -> edist(x, year))
    res = optimize(f,-50,50)
    days_per[i] = Optim.minimizer(res)[1]
end

plot((years), days_eq, legend=false, dpi=150)
xlabel!("Year")
ylabel!("Day in March")
title!("Date of vernal equinox")
savefig("equinox_dates.png")

plot((years), days_per, legend=false, dpi=150)
xlabel!("Year")
ylabel!("Day in Jan")
title!("Date of perihelion")
savefig("perihelion_dates.png")

years = -100000:100:100000 #hide
days_eq = zeros(length(years)) #hide
for (i,year) in enumerate(years) #hide
    f = (x -> zdiff(x, year)) #hide
    days_eq[i] = find_zeros(f,-30,60)[1] #hide
end #hide

plot((years / 1000), days_eq, legend=false, dpi=150) #hide
xlabel!("kyr") #hide
ylabel!("Day in March") #hide
title!("Date of vernal equinox") #hide
savefig("equinox_dates_long.png") #hide
```
![](equinox_dates.png)
![](perihelion_dates.png)

## Variations in date of vernal equinox on millenial timescales
![](equinox_dates_long.png)

# Insolation Examples

## Diurnal cycle of insolation
```@example
using Dates

include("plot_diurnal_cycle.jl")

# Pasadena in January
lat, lon = [34.15, -118.14]
date = DateTime(2020, 01, 10)
timezone = +8 # Pacific Standard Time
diurnal_cycle(lat, lon, date, timezone, "Pasadena_January.png")

# Finland in June
lat, lon = [66.50, 25.73]
date = DateTime(2020, 06, 10)
timezone = -3 # Eastern European Summer Time
diurnal_cycle(lat, lon, date, timezone, "Finland_June.png")
```
![](Pasadena_January.png)
![](Finland_June.png)

## Daily-mean insolation
### Insolation in J2000
```@example
using CLIMAParameters
using CLIMAParameters.Planet
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

include("plot_insolation.jl")

γ0 = obliq_epoch(param_set)
ϖ0 = lon_perihelion_epoch(param_set)
e0 = eccentricity_epoch(param_set)

days, lats, F0 = calc_day_lat_insolation(365, 180, param_set)
title = format("γ = {:.2f}°, ϖ = {:.2f}°, e = {:.2f}", rad2deg(γ0), rad2deg(ϖ0), e0) #hide
plot_day_lat_insolation(days, lats, F0, "YlOrRd", title, "insol_example1.png")
```
![](insol_example1.png)

### Insolation with smaller obliquity
```@example
using CLIMAParameters # hide
using CLIMAParameters.Planet # hide
struct EarthParameterSet <: AbstractEarthParameterSet end # hide
const param_set = EarthParameterSet() # hide
include("plot_insolation.jl") # hide
γ0 = obliq_epoch(param_set) # hide
ϖ0 = lon_perihelion_epoch(param_set) # hide
e0 = eccentricity_epoch(param_set) # hide
days, lats, F0 = calc_day_lat_insolation(365, 180, param_set) # hide

# decrease γ to 20.0°
CLIMAParameters.Planet.obliq_epoch(::EarthParameterSet) = deg2rad(20.0)
γ1 = obliq_epoch(param_set)
days, lats, F2 = calc_day_lat_insolation(365, 180, param_set)

title = format("γ = {:.2f}°, ϖ = {:.2f}°, e = {:.2f}", rad2deg(γ1), rad2deg(ϖ0), e0) # hide
plot_day_lat_insolation(days,lats,F2,"YlOrRd",  title, "insol_example2a.png")
title = format("insolation diff: γ' = {:.2f}° - γ = {:.2f}°", rad2deg(γ1), rad2deg(γ0)) # hide
plot_day_lat_insolation(days, lats, F2-F0, "PRGn", title, "insol_example2b.png")
```
![](insol_example2a.png)
![](insol_example2b.png)

### Insolation with very large obliquity (like Uranus)
```@example
using CLIMAParameters # hide
using CLIMAParameters.Planet # hide
struct EarthParameterSet <: AbstractEarthParameterSet end # hide
const param_set = EarthParameterSet() # hide
include("plot_insolation.jl") # hide
γ0 = obliq_epoch(param_set) # hide
ϖ0 = lon_perihelion_epoch(param_set) # hide
e0 = eccentricity_epoch(param_set) # hide
days, lats, F0 = calc_day_lat_insolation(365, 180, param_set) # hide

# now change obliquity to 97.86°
CLIMAParameters.Planet.obliq_epoch(::EarthParameterSet) = deg2rad(97.86)
γ4 = obliq_epoch(param_set)
days, lats, F5 = calc_day_lat_insolation(365, 180, param_set)

title = format("γ = {:.2f}°, ϖ = {:.2f}°, e = {:.2f}", rad2deg(γ4), rad2deg(ϖ0), e0) # hide
plot_day_lat_insolation(days,lats,F5,"YlOrRd", title, "insol_example3a.png")
title = format("insolation diff: γ' = {:.2f}° - γ = {:.2f}°", rad2deg(γ4), rad2deg(γ0)) # hide
plot_day_lat_insolation(days, lats, F5-F0, "PRGn", title, "insol_example3b.png")

```
![](insol_example3a.png)
![](insol_example3b.png)
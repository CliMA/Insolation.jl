# Insolation Examples

This page demonstrates various use cases of `Insolation.jl` through visualizations and practical examples. These examples show how to calculate and visualize solar radiation patterns across different timescales and locations.

## Diurnal Cycle of Insolation

The diurnal cycle shows how insolation and solar zenith angle vary throughout a day at a specific location. This is useful for understanding the daily energy input and for applications such as solar power forecasting.

**Pasadena in Winter**: Clear diurnal cycle with peak insolation at solar noon. Solar zenith angle reaches minimum (~35°) at midday.

```@example diurnal
using Insolation
using Dates

include("plot_diurnal_cycle.jl")

# Load orbital data
od = OrbitalDataSplines()

# Example 1: Pasadena, California in January (mid-latitude winter)
lat, lon = [34.15, -118.14]  # Pasadena coordinates
date = DateTime(2020, 01, 10)
timezone = +8  # Pacific Standard Time (UTC+8)
diurnal_cycle(lat, lon, date, od, timezone, "Pasadena_January.png")
nothing # hide
```
![](Pasadena_January.png)

**Finland in Summer**: At the Arctic Circle, the sun does not set at northern summer solitice, resulting in continuous daylight and non-zero insolation for the 24-hour period. This "midnight sun" phenomenon is visible as the sustained insolation during night hours.

```@example diurnal
# Example 2: Rovaniemi, Finland in June (Arctic summer - midnight sun)
lat, lon = [66.50, 25.73]  # Rovaniemi coordinates (Arctic Circle)
date = DateTime(2020, 06, 20)  # Summer solstice
timezone = -3  # Eastern European Summer Time (UTC+3)
diurnal_cycle(lat, lon, date, od, timezone, "Finland_June.png")
nothing # hide
```
![](Finland_June.png)

## Latitudinal and Seasonal Variations

The following examples show how daily-averaged insolation varies with latitude and day of year, revealing Earth's seasonal cycles and the role of orbital parameters.

### Modern Climate (J2000 Epoch)
```@example insolation_examples
import Insolation
import Insolation.Parameters as IP
import ClimaParams as CP

FT = Float64
param_set = IP.InsolationParameters(FT)

include("plot_insolation.jl")

# Get current epoch orbital parameters
γ0 = IP.obliq_epoch(param_set)
ϖ0 = IP.lon_perihelion_epoch(param_set)
e0 = IP.eccentricity_epoch(param_set)
od = Insolation.OrbitalDataSplines()

# Calculate daily-mean insolation across latitudes and days
days, lats, F0 = calc_day_lat_insolation(od, 365, 180, param_set)
title = format("Modern Earth: γ = {:.2f}°, ϖ = {:.2f}°, e = {:.3f}", rad2deg(γ0), rad2deg(ϖ0), e0) #hide
plot_day_lat_insolation(days, lats, F0, "YlOrRd", title, "insol_example1.png")
nothing # hide
```

This plot shows daily-averaged TOA insolation as a function of latitude and day of year for modern Earth (J2000 epoch). Key features:
- **Solstices**: Maximum insolation at high latitudes during local summer (June in NH, December in SH)
- **Equinoxes**: Symmetric insolation distribution around March 21 and September 21
- **Annual Mean**: Right panel shows the latitudinal gradient, with maximum near the equator

![](insol_example1.png)

### Effect of Reduced Obliquity

Obliquity (axial tilt) controls the strength of seasonal cycles. Here we reduce obliquity from 23.44° to 20° to demonstrate its impact.
```@example insolation_examples
# Reduce obliquity to 20.0° (from current 23.44°)
param_set_low_obliq = IP.InsolationParameters(FT, (; obliq_epoch = deg2rad(20.0)))
γ1 = IP.obliq_epoch(param_set_low_obliq)
days, lats, F2 = calc_day_lat_insolation(od, 365, 180, param_set_low_obliq)

title = format("Low Obliquity: γ = {:.2f}°, ϖ = {:.2f}°, e = {:.3f}", rad2deg(γ1), rad2deg(ϖ0), e0) # hide
plot_day_lat_insolation(days,lats,F2,"YlOrRd",  title, "insol_example2a.png")
title = format("Difference: γ = {:.2f}° minus γ = {:.2f}°", rad2deg(γ1), rad2deg(γ0)) # hide
plot_day_lat_insolation(days, lats, F2-F0, "PRGn", title, "insol_example2b.png")
nothing # hide
```

**Absolute Insolation** (top): With reduced obliquity, polar regions receive less summer insolation, and the tropical regions receive slightly more year-round.

![](insol_example2a.png)

**Difference Plot** (bottom): Shows that reducing obliquity:
- **Decreases** high-latitude summer insolation (blue regions)
- **Increases** tropical insolation slightly (purple regions)
- Results in **weaker seasonal cycles** overall

This demonstrates that lower obliquity leads to milder seasons - important for understanding climates on other planets or Earth's past climates.

![](insol_example2b.png)

### Extreme Obliquity (Uranus-like Configuration)

Uranus has an extreme axial tilt of 97.86°, essentially rotating "on its side." This creates dramatically different seasonal patterns.
```@example insolation_examples
# Set obliquity to 97.86° (Uranus's axial tilt)
param_set_uranus = IP.InsolationParameters(FT, (;obliq_epoch = deg2rad(97.86)))
γ4 = IP.obliq_epoch(param_set_uranus)
days, lats, F5 = calc_day_lat_insolation(od, 365, 180, param_set_uranus)

title = format("Uranus Obliq.: γ = {:.2f}°, ϖ = {:.2f}°, e = {:.3f}", rad2deg(γ4), rad2deg(ϖ0), e0) # hide
plot_day_lat_insolation(days,lats,F5,"YlOrRd", title, "insol_example3a.png")
title = format("Difference: γ = {:.2f}° minus γ = {:.2f}°", rad2deg(γ4), rad2deg(γ0)) # hide
plot_day_lat_insolation(days, lats, F5-F0, "PRGn", title, "insol_example3b.png")
nothing # hide
```

**Absolute Insolation** (top): With extreme obliquity, the insolation pattern is dramatically different:
- **Polar regions** receive intense summer insolation, exceeding tropical values
- **Equatorial regions** have two insolation peaks per year
- **Extreme seasons** with long polar day/night periods

![](insol_example3a.png)

**Difference Plot** (bottom): Compared to Earth's current configuration:
- Polar summers receive **much more** insolation (>500 W/m² increase)
- Annual-mean distribution shifts toward maxima at the poles

This extreme case illustrates how orbital parameters fundamentally shape a planet's climate. Such configurations help us understand exotic exoplanets and the range of possible planetary climates.

![](insol_example3b.png)

## Summary

These examples demonstrate that `Insolation.jl` can:
- Calculate diurnal cycles for any location and date
- Visualize latitudinal and seasonal patterns
- Explore sensitivity to orbital parameters
- Model both Earth-like and exotic planetary configurations

For paleoclimate applications with time-varying orbital parameters, see [Milankovitch Cycles](@ref).
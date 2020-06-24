# Solar Zenith Angle Examples

## Example 1
```@example
include("plot_zenith.jl")

# 2018 orbital constants
γ0 = 23.44
ϖ0 = 282.95
e0 = 0.017

days, lats, S = calc_day_lat_sza(365, 180, γ0, ϖ0, e0)
title = "obliq=" * "$(γ0)" * ", perihelion=" * "$(ϖ0)" * ", ecc=" * "$(e0)"
plot_day_lat_sza(days, lats, S, "YlOrRd", title, "sza_example1.png")
```
![](sza_example1.png)

## Example 2
```@example
include("plot_zenith.jl") # hide

# turn longitude of perihelion by 180°
γ0 = 23.44
ϖ0 = 282.95 # hide
ϖ1 = ϖ0 + 180.0
e0 = 0.017

days, lats, S0 = calc_day_lat_sza(365, 180, γ0, ϖ0, e0) # hide
days, lats, S1 = calc_day_lat_sza(365, 180, γ0, ϖ1, e0)

title = "obliq=" * "$(γ0)" * ", perihelion=" * "$(ϖ1)" * ", ecc=" * "$(e0)"
plot_day_lat_sza(days, lats, S1, "YlOrRd",  title, "sza_example2a.png")

title = "SZA diff: perihelion0=" * "$(ϖ0)" * ", perihelion1=" * "$(ϖ1)"
plot_day_lat_sza(days,lats,S1-S0,"PRGn", title, "sza_example2b.png")
```
![](sza_example2a.png)
![](sza_example2b.png)

## Example 3
```@example
include("plot_zenith.jl") # hide

# perihelion back to normal. decrease γ to 22.0°
γ0 = 23.44 # hide
γ1 = 22.0
ϖ0 = 282.95
e0 = 0.017

days, lats, S0 = calc_day_lat_sza(365, 180, γ0, ϖ0, e0) # hide
days, lats, S2 = calc_day_lat_sza(365, 180, γ1, ϖ0, e0)

title = "obliq=" * "$(γ1)" * ", perihelion=" * "$(ϖ0)" * ", ecc=" * "$(e0)"
plot_day_lat_sza(days,lats,S2,"YlOrRd",  title, "sza_example3a.png")

title = "SZA diff: obliq0=" * "$(γ0)" * ", obliq1=" * "$(γ1)"
plot_day_lat_sza(days,lats,S2-S0,"PRGn", title, "sza_example3b.png")
```
![](sza_example3a.png)
![](sza_example3b.png)

## Example 4
```@example
include("plot_zenith.jl") # hide

# decrease γ further to 18.0°
γ0 = 23.44 # hide
γ2 = 18.0
ϖ0 = 282.95
e0 = 0.017

days, lats, S0 = calc_day_lat_sza(365, 180, γ0, ϖ0, e0) # hide
days, lats, S3 = calc_day_lat_sza(365, 180, γ2, ϖ0, e0)

title = "obliq=" * "$(γ2)" * ", perihelion=" * "$(ϖ0)" * ", ecc=" * "$(e0)"
plot_day_lat_sza(days,lats,S3,"YlOrRd", title, "sza_example4a.png")

title = "SZA diff: obliq0=" * "$(γ0)" * ", obliq1=" * "$(γ2)"
plot_day_lat_sza(days,lats,S3-S0,"PRGn", title, "sza_example4b.png")
```
![](sza_example4a.png)
![](sza_example4b.png)

## Example 5
```@example
include("plot_zenith.jl") # hide

# now change obliquity to 60.0°
γ0 = 23.44 # hide
γ3 = 60.0
ϖ0 = 282.95
e0 = 0.017

days, lats, S0 = calc_day_lat_sza(365, 180, γ0, ϖ0, e0) # hide
days, lats, S4 = calc_day_lat_sza(365, 180, γ3, ϖ0, e0)

title = "obliq=" * "$(γ3)" * ", perihelion=" * "$(ϖ0)" * ", ecc=" * "$(e0)"
plot_day_lat_sza(days,lats,S4,"YlOrRd", title, "sza_example5a.png")

title = "SZA diff: obliq0=" * "$(γ0)" * ", obliq1=" * "$(γ3)"
plot_day_lat_sza(days,lats,S4-S0,"PRGn", title, "sza_example5b.png")
```
![](sza_example5a.png)
![](sza_example5b.png)

## Example 6
```@example
include("plot_zenith.jl") # hide

# now change obliquity to 97.86°
γ0 = 23.44 # hide
γ4 = 97.86
ϖ0 = 282.95
e0 = 0.017

days, lats, S0 = calc_day_lat_sza(365, 180, γ0, ϖ0, e0) # hide
days, lats, S5 = calc_day_lat_sza(365, 180, γ4, ϖ0, e0)

title = "obliq=" * "$(γ4)" * ", perihelion=" * "$(ϖ0)" * ", ecc=" * "$(e0)"
plot_day_lat_sza(days,lats,S5,"YlOrRd", title, "sza_example6a.png")

title = "SZA diff: obliq0=" * "$(γ0)" * ", obliq1=" * "$(γ4)"
plot_day_lat_sza(days,lats,S5-S0,"PRGn", title, "sza_example6b.png")

```
![](sza_example6a.png)
![](sza_example6b.png)
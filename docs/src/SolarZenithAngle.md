# Solar Zenith Angle

```@meta
CurrentModule = Insolation.SolarZenithAngle
```

## Example 1
```@example
include("../../examples/plot_insolation.jl")

# 2018 orbital constants
γ0 = 23.44
ϖ0 = 282.95
e0 = 0.017

days, lats, F0 = calc_day_lat_insolation(365, 180, γ0, ϖ0, e0)
title = "obliq=" * "$(γ0)" * ", perihelion=" * "$(ϖ0)" * ", ecc=" * "$(e0)"
plot_day_lat_insolation(days, lats, F0, "YlOrRd", title, "example1.png")
```
![](example1.png)

## Example 2
```@example
include("../../examples/plot_insolation.jl") # hide

# turn longitude of perihelion by 180°
γ0 = 23.44
ϖ0 = 282.95 # hide
ϖ1 = ϖ0 + 180.0
e0 = 0.017

days, lats, F0 = calc_day_lat_insolation(365, 180, γ0, ϖ0, e0) # hide
days, lats, F1 = calc_day_lat_insolation(365, 180, γ0, ϖ1, e0)

title = "obliq=" * "$(γ0)" * ", perihelion=" * "$(ϖ1)" * ", ecc=" * "$(e0)"
plot_day_lat_insolation(days, lats, F1, "YlOrRd",  title, "example2a.png")

title = "insolation diff: perihelion0=" * "$(ϖ0)" * ", perihelion1=" * "$(ϖ1)"
plot_day_lat_insolation(days,lats,F1-F0,"PRGn", title, "example2b.png")
```
![](example2a.png)
![](example2b.png)

## Example 3
```@example
include("../../examples/plot_insolation.jl") # hide

# perihelion back to normal. decrease γ to 22.0°
γ0 = 23.44 # hide
γ1 = 22.0
ϖ0 = 282.95
e0 = 0.017

days, lats, F0 = calc_day_lat_insolation(365, 180, γ0, ϖ0, e0) # hide
days, lats, F2 = calc_day_lat_insolation(365, 180, γ1, ϖ0, e0)

title = "obliq=" * "$(γ1)" * ", perihelion=" * "$(ϖ0)" * ", ecc=" * "$(e0)"
plot_day_lat_insolation(days,lats,F2,"YlOrRd",  title, "example3a.png")

title = "insolation diff: obliq0=" * "$(γ0)" * ", obliq1=" * "$(γ1)"
plot_day_lat_insolation(days,lats,F2-F0,"PRGn", title, "example3b.png")
```
![](example3a.png)
![](example3b.png)

## Example 4
```@example
include("../../examples/plot_insolation.jl") # hide

# decrease γ further to 18.0°
γ0 = 23.44 # hide
γ2 = 18.0
ϖ0 = 282.95
e0 = 0.017

days, lats, F0 = calc_day_lat_insolation(365, 180, γ0, ϖ0, e0) # hide
days, lats, F3 = calc_day_lat_insolation(365, 180, γ2, ϖ0, e0)

title = "obliq=" * "$(γ2)" * ", perihelion=" * "$(ϖ0)" * ", ecc=" * "$(e0)"
plot_day_lat_insolation(days,lats,F3,"YlOrRd", title, "example4a.png")

title = "insolation diff: obliq0=" * "$(γ0)" * ", obliq1=" * "$(γ2)"
plot_day_lat_insolation(days,lats,F3-F0,"PRGn", title, "example4b.png")
```
![](example4a.png)
![](example4b.png)

## Example 5
```@example
include("../../examples/plot_insolation.jl") # hide

# now change obliquity to 60.0°
γ0 = 23.44 # hide
γ3 = 60.0
ϖ0 = 282.95
e0 = 0.017

days, lats, F0 = calc_day_lat_insolation(365, 180, γ0, ϖ0, e0) # hide
days, lats, F4 = calc_day_lat_insolation(365, 180, γ3, ϖ0, e0)

title = "obliq=" * "$(γ3)" * ", perihelion=" * "$(ϖ0)" * ", ecc=" * "$(e0)"
plot_day_lat_insolation(days,lats,F4,"YlOrRd", title, "example5a.png")

title = "insolation diff: obliq0=" * "$(γ0)" * ", obliq1=" * "$(γ3)"
plot_day_lat_insolation(days,lats,F4-F0,"PRGn", title, "example5b.png")
```
![](example5a.png)
![](example5b.png)

## Example 6
```@example
include("../../examples/plot_insolation.jl") # hide

# now change obliquity to 97.86°
γ0 = 23.44 # hide
γ4 = 97.86
ϖ0 = 282.95
e0 = 0.017

days, lats, F0 = calc_day_lat_insolation(365, 180, γ0, ϖ0, e0) # hide
days, lats, F5 = calc_day_lat_insolation(365, 180, γ4, ϖ0, e0)

title = "obliq=" * "$(γ4)" * ", perihelion=" * "$(ϖ0)" * ", ecc=" * "$(e0)"
plot_day_lat_insolation(days,lats,F5,"YlOrRd", title, "example6a.png")

title = "insolation diff: obliq0=" * "$(γ0)" * ", obliq1=" * "$(γ4)"
plot_day_lat_insolation(days,lats,F5-F0,"PRGn", title, "example6b.png")

```
![](example6a.png)
![](example6b.png)
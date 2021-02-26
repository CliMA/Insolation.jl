# Insolation Examples
These examples are from a homework assignment in ESE 101 at Caltech

## Example 1
```@example
include("plot_insolation.jl")

# orbital constants
γ0 = γ_epoch()
ϖ0 = ϖ_epoch()
e0 = e_epoch()

days, lats, F0 = calc_day_lat_insolation(365, 180, γ0, ϖ0, e0)
title = format("g = {:.2f}, w = {:.2f}, e = {:.2f}", γ0, ϖ0, e0) #hide
plot_day_lat_insolation(days, lats, F0, "YlOrRd", title, "insol_example1.png")
```
![](insol_example1.png)

## Example 2
```@example
include("plot_insolation.jl") # hide

# turn longitude of perihelion by 180°
γ0 = γ_epoch()
ϖ0 = ϖ_epoch() # hide
ϖ1 = ϖ0 + π
e0 = e_epoch()

days, lats, F0 = calc_day_lat_insolation(365, 180, γ0, ϖ0, e0) # hide
days, lats, F1 = calc_day_lat_insolation(365, 180, γ0, ϖ1, e0)

title = format("g = {:.2f}, w = {:.2f}, e = {:.2f}", γ0, ϖ1, e0) # hide
plot_day_lat_insolation(days, lats, F1, "YlOrRd",  title, "insol_example2a.png")
title = format("insolation diff: w0 = {:.2f}, w1 = {:.2f}", ϖ0, ϖ1) # hide
plot_day_lat_insolation(days, lats, F1-F0, "PRGn", title, "insol_example2b.png")
```
![](insol_example2a.png)
![](insol_example2b.png)

## Example 3
```@example
include("plot_insolation.jl") # hide

# perihelion back to normal. decrease γ to 22.0°
γ0 = γ_epoch() # hide
γ1 = deg2rad(22.0)
ϖ0 = ϖ_epoch()
e0 = e_epoch()

days, lats, F0 = calc_day_lat_insolation(365, 180, γ0, ϖ0, e0) # hide
days, lats, F2 = calc_day_lat_insolation(365, 180, γ1, ϖ0, e0)

title = format("g = {:.2f}, w = {:.2f}, e = {:.2f}", γ1, ϖ0, e0) # hide
plot_day_lat_insolation(days,lats,F2,"YlOrRd",  title, "insol_example3a.png")
title = format("insolation diff: g0 = {:.2f}, g1 = {:.2f}", γ0, γ1) # hide
plot_day_lat_insolation(days, lats, F2-F0, "PRGn", title, "insol_example3b.png")
```
![](insol_example3a.png)
![](insol_example3b.png)

## Example 4
```@example
include("plot_insolation.jl") # hide

# decrease γ further to 18.0°
γ0 = γ_epoch() # hide
γ2 = deg2rad(18.0)
ϖ0 = ϖ_epoch()
e0 = e_epoch()

days, lats, F0 = calc_day_lat_insolation(365, 180, γ0, ϖ0, e0) # hide
days, lats, F3 = calc_day_lat_insolation(365, 180, γ2, ϖ0, e0)

title = format("g = {:.2f}, w = {:.2f}, e = {:.2f}", γ2, ϖ0, e0) # hide
plot_day_lat_insolation(days,lats,F3,"YlOrRd", title, "insol_example4a.png")
title = format("insolation diff: g0 = {:.2f}, g1 = {:.2f}", γ0, γ2) # hide
plot_day_lat_insolation(days, lats, F3-F0, "PRGn", title, "insol_example4b.png")
```
![](insol_example4a.png)
![](insol_example4b.png)

## Example 5
```@example
include("plot_insolation.jl") # hide

# now change obliquity to 60.0°
γ0 = γ_epoch() # hide
γ3 = deg2rad(60.0)
ϖ0 = ϖ_epoch()
e0 = e_epoch()

days, lats, F0 = calc_day_lat_insolation(365, 180, γ0, ϖ0, e0) # hide
days, lats, F4 = calc_day_lat_insolation(365, 180, γ3, ϖ0, e0)

title = format("g = {:.2f}, w = {:.2f}, e = {:.2f}", γ3, ϖ0, e0) # hide
plot_day_lat_insolation(days,lats,F4,"YlOrRd", title, "insol_example5a.png")
title = format("insolation diff: g0 = {:.2f}, g1 = {:.2f}", γ0, γ3) # hide
plot_day_lat_insolation(days, lats, F4-F0, "PRGn", title, "insol_example5b.png")
```
![](insol_example5a.png)
![](insol_example5b.png)

## Example 6
```@example
include("plot_insolation.jl") # hide

# now change obliquity to 97.86°
γ0 = γ_epoch() # hide
γ4 = deg2rad(97.86)
ϖ0 = ϖ_epoch()
e0 = e_epoch()

days, lats, F0 = calc_day_lat_insolation(365, 180, γ0, ϖ0, e0) # hide
days, lats, F5 = calc_day_lat_insolation(365, 180, γ4, ϖ0, e0)

title = format("g = {:.2f}, w = {:.2f}, e = {:.2f}", γ4, ϖ0, e0) # hide
plot_day_lat_insolation(days,lats,F5,"YlOrRd", title, "insol_example6a.png")
title = format("insolation diff: g0 = {:.2f}, g1 = {:.2f}", γ0, γ4) # hide
plot_day_lat_insolation(days, lats, F5-F0, "PRGn", title, "insol_example6b.png")

```
![](insol_example6a.png)
![](insol_example6b.png)
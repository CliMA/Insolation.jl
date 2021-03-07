# Insolation.jl

```@meta
CurrentModule = Insolation
```

`Insolation.jl` is a library that calculates the zenith angle and insolation 
at a given point on Earth (lat, lon) for a given date/time and orbital configuration 
(obliquity, eccentricity, and longitude of perihelion). The library is split
between two files, `ZenithAngleCalc.jl` which calculates the zenith angle, azimuth angle, and Earth-sun distance, 
and `InsolationCalc.jl` which calculates the insolation given a zenith angle.

The zenith angle and insolation can both be calculated either as instantaneous 
values or as daily averaged values. The functions in `ZenithAngleCalc.jl` are 
overwritten to accept a variety of inputs, either calculating the orbital parameters 
given a DateTime object, or prescribing the orbital parameters as input.

The equations used for calculating orbital parameters are from Tapio Schneider's textbook draft. 
See [Zenith Angle Equations](@ref) for more details.

This package has been designed to be flexible for non-Earth settings as well.
This specific functionality has not yet been implemented [as of 22 Feb 2021].

## Authors
`Insolation.jl` is being developed by [the Climate Modeling Alliance](https://clima.caltech.edu).
Specifically it has been developed to be used by [`RRTMGP.jl`](https://github.com/CliMA/RRTMGP.jl) 
and the [`ClimateMachine.jl`](https://github.com/CliMA/ClimateMachine.jl).
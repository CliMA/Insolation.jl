# Zenith Angle Equations
These equations are from Tapio's textbook draft chapter 3.

## Mean Anomaly at Vernal Equinox (3.11)
```math
\begin{aligned}
\beta &= (1-e^2)^{1/2} \\
M_{VE} &= - \varpi + (e + e^{3/4})(1+\beta) \sin(\varpi) - \frac{e^2}{2} \left( \frac{1}{2} + \beta \right) \sin(2 \varpi) + \frac{e^3}{4} \left( \frac{1}{3} + \beta \right) \sin(3 \varpi)
\end{aligned}
```

## Mean Anomaly (3.10)
```math
M = \frac{2\pi (t-t_{VE})}{Y_a} + M_{VE}
```
where $t_{VE}$ is the time of the vernal equinox.

## True Anomaly (3.8)
```math
A = M + (2e - e^{3/4}) \sin(M) + \frac{5}{4} e^2 \sin(2M) + \frac{13}{12} e^3 \sin(3M)
```

## True Longitude (3.9)
```math
L = A + \varpi
```

## Declination
```math
\sin(\delta) = \sin(\gamma) \sin(L)
```

## Hour Angle
```math
\eta = \frac{2\pi (t-t_s)}{T_d}
```
where $t_s$ is the time of local solar noon and $T_d$ is the length of a day.

## Zenith Angle
```math
\cos(\theta) = \cos(\phi) \cos(\delta) \cos(\eta) + \sin(\phi) \sin(\delta)
```

## Azimuth Angle
```math
\tan(\zeta) = \frac{\sin(\eta)}{\cos(\eta)\sin(\phi) - \tan(\delta)\cos(\phi)}
```

## Earth-Sun Distance
```math
d = \frac{1-e^2}{1+e\cos(A)} d_0
```
where $d_0$ is 1 AU.
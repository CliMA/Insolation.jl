# Zenith Angle Equations
These equations are from the book "Astronomical Algorithms"
by Jean Meeus.

## Times: Julian Date/Julian Century
The Julian century (``jc``) is defined as the centuries
since 1 Jan 2000 12:00 UTC.

## True Solar Longitude/Anomaly
```math
\begin{aligned}
ML &= 280.46646 + 36000.76983 jc + 0.0003032 jc^2 \\
MA &= 357.52911 + 35999.05029 jc - 0.0001537 jc^2 \\
SC &= \sin(MA) (1.914602-0.004817 jc-0.000014 jc^2) \\
&\hspace{1cm}+ \sin(2 MA) (0.019993-0.000101 jc) \\
&\hspace{1cm}+ \sin(3 MA) (0.000289) \\
TL &= ML + SC \\
TA &= MA + SC
\end{aligned}
```

## Eccentricity and Obliquity
```math
\begin{aligned}
e &= 0.016708634 - 0.000042037 jc - 0.0000001267 jc^2 \\
\gamma &= 23.439291 - 0.01300417 jc - 1.638889e-7 jc^2 + 5.036111e-7 jc^3
\end{aligned}
```

## Declination and Right Acension
```math
\begin{aligned}
\delta &= \arcsin(sin(\gamma) \sin(TL)) \\
RA &= \arctan \left( \frac{\cos(\gamma) \sin(TL)}{\cos(TL)} \right)
\end{aligned}
```

## Greenwich mean sidereal time
```math
GMST = 6.6974243242 + 2400.117188 jc + UTC_{hours}
```

## Hour Angle
```math
\eta = GMST + \lambda - RA
```

## Zenith Angle
```math
\theta = \arccos(\cos(\phi) \cos(\delta) \cos(\eta) + \sin(\phi) \sin(\delta))
```

## Daily averaged zenith angle
```math
\begin{aligned}
T &= \tan(\phi) \tan(\delta) \\
\eta_d &= \left\{ \begin{array}{ll} \eta_d = \pi & T \geq 1 \\ \eta_d = 0 & T \leq -1 \\ \eta_d = \arccos(-T) & else \end{array} \right. \\
\overline{\theta} &= \arccos \left( \left(\frac{1}{\pi} \right) (\eta_d \sin(\phi) \sin(\delta) + \cos(\phi) \cos(\delta) \sin(\eta_d)) \right)
\end{aligned}
```

## Earth-Sun Distance
```math
d_{au} = \frac{1.000001018 (1.0 - e^2)}{1.0 + e \cos(TA)}
```
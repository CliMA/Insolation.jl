# Zenith Angle Equations
These equations are from Tapio Schneiders's textbook draft chapter 3.

## Mean Anomaly (3.6)
```math
\begin{aligned}
\beta &= \sqrt{1-e^2} \\
M_{V0} &= -\varpi_0 + \left(e+\frac{1}{4}e^3 \right) (1+\beta)sin(\varpi_0) - \frac{1}{2}e^2 \left(\frac{1}{2}+\beta \right) sin(2\varpi_0) + \frac{1}{4}e^3 \left(\frac{1}{3}+\beta \right) sin(3\varpi_0) \\
t_V &= \frac{Y_a (M_{V0} - M_0)}{2\pi} + t_0 \\
M_V &= -\varpi + \left(e+\frac{1}{4}e^3 \right) (1+\beta)sin(\varpi) - \frac{1}{2}e^2 \left(\frac{1}{2}+\beta \right) sin(2\varpi) + \frac{1}{4}e^3 \left(\frac{1}{3}+\beta  \right) sin(3\varpi) \\
M &= \frac{2\pi (t - t_V)}{Y_a} + M_V
\end{aligned}
```
where $e$ is the oribital eccentricity, $t_0$ is the time at the epoch (J2000), defined as January 1, 2000 at 12hr UTC, 
$M_0$ is the mean anomaly at the epoch, $\varpi_0$ is the longitude of perihelion at the epoch,
and $Y_a$ is the length of the anomalistic year.

## True Anomaly (3.8)
```math
A = M + (2e - e^{3/4}) \sin(M) + \frac{5}{4} e^2 \sin(2M) + \frac{13}{12} e^3 \sin(3M)
```

## True Longitude (3.9)
```math
L = A + \varpi
```

## Declination (3.16)
```math
\sin \delta = \sin \gamma \sin L
```
where $\gamma$ is the orbital obliquity.

## Hour Angle (3.17)
```math
\eta = \frac{2\pi (t-t_s)}{T_d}
```
where $t_s$ is the time of local solar noon and $T_d$ is the length of a day.

## Zenith Angle (3.18)
```math
\cos \theta = \cos \phi \cos \delta \cos \eta + \sin \phi \sin \delta
```
where $\phi$ is the latitude.

## Sunrise/Sunset Angle (3.19)
```math
\cos \eta_d = - \tan \phi \tan \delta
```

## Daily-averaged Zenith Angle (3.20)
```math
\overline{\cos \theta} = \frac{1}{\pi} \left( \eta_d \sin \phi \sin \delta + \cos \phi \cos \delta \cos \eta_d \right)
```

## Azimuth Angle
```math
\tan \zeta = \frac{\sin \eta}{\cos \eta \sin \phi - \tan \delta \cos \phi}
```

## Earth-Sun Distance (3.1)
```math
d = \frac{1-e^2}{1+e\cos A} d_0
```
where $d_0$ is 1 AU.
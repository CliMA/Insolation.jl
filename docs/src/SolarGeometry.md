# Solar Geometry and Insolation

This page provides the mathematical formulations used in `Insolation.jl` to calculate solar geometry and insolation. The equations are from Tapio Schneider and Lenka Novak's textbook draft ["Physics of Earth's Climate"](https://climate-dynamics.org/wp-content/uploads/2017/04/Climate_Book.pdf), Chapter 3.

## Overview

Calculating solar position in the sky and insolation requires a sequence of astronomical computations that determine:

1. The planet's position in its orbit (mean and true anomaly)
2. The latitude of the subsolar point (declination)
3. The longitude relative to solar noon (hour angle)
4. The geometric relationships between the sun and a surface location (zenith and azimuth angles)
5. The planet-star distance for radiation calculations

## Mean Anomaly

The **mean anomaly** $M$ represents the angular position of a planet in its orbit, assuming uniform circular motion. For an elliptical orbit, this is a convenient starting point for more accurate calculations.

The mean anomaly at current time $t$ is

```math
M = \frac{2\pi (t - t_0)}{Y_a} + M_0,
```

where we have:

- The time at the epoch (J2000), $t_0$, typically defined as January 1, 2000 at 12:00 Terrestrial Time (TT), corresponding to 11:59 UTC
- The mean anomaly at the epoch, $M_0$ [radians]
- The length of the anomalistic year, $Y_a$ (period from perihelion to perihelion) [seconds]

The mean anomaly increases linearly with time at a rate of $2\pi/Y_a$ radians per second.

## True Anomaly

The **true anomaly** $A$ is the actual angular position of the planet in its elliptical orbit, measured from perihelion (the point of closest approach to the star). Unlike the mean anomaly, the true anomaly accounts for the planet's varying orbital speed due to Kepler's laws.

The true anomaly is computed from the mean anomaly using a series expansion,

```math
A = M + \left( 2e - \frac{1}{4}e^{3} \right) \sin(M) + \frac{5}{4} e^2 \sin(2M) + \frac{13}{12} e^3 \sin(3M) + \mathcal{O}(e^4),
```

where $e$ is the orbital eccentricity. This series approximation is accurate to order $e^3$ and is sufficient for Earth's relatively circular orbit ($e \approx 0.017$).

!!! warning "Low-eccentricity approximation"
    The series error grows rapidly with eccentricity (exceeding a few degrees near $e \approx 0.5$). For highly eccentric orbits, it should be replaced by an exact solution of Kepler's equation (e.g., a Newton–Raphson iteration for the eccentric anomaly). Because the true anomaly feeds the declination, planet-star distance, and equation of time, all of these inherit the same eccentricity limitation.

## Solar Longitude

The **solar longitude** $L_s$ (also called ecliptic longitude or true longitude) is the ecliptic longitude of the Sun, measured from the vernal equinox. It combines the orbital phase (true anomaly $A$) with the longitude of perihelion $\varpi$: the Sun's ecliptic longitude at perihelion, measured from the vernal equinox (≈ 283° for present-day Earth):

```math
L_s = A + \varpi.
```

For Earth, $\varpi$ varies slowly due to precession (period $\sim 21{,}000$ years).

## Declination

The **solar declination** $\delta$ is the angle between the sun's rays and the equatorial plane, or the latitude of the subsolar point (where the sun is in zenith at solar noon). It determines the subsolar latitude (where the sun is directly overhead at solar noon). The declination varies between $\pm\gamma$ (obliquity) over the course of a year.

The sine of the declination angle is

```math
\sin \delta = \sin \gamma \sin L_s,
```

where we have:

- The orbital obliquity $\gamma$ (axial tilt) [radians]
- The solar longitude $L_s$ [radians]

For Earth's current obliquity of $\gamma \approx 23.44°$, the declination ranges from $-23.44°$ (winter solstice) to $+23.44°$ (summer solstice), passing through $0°$ at the equinoxes.

## Equation of Time

The **equation of time** corrects for the difference between **apparent solar time** (when the sun is actually at its highest point) and **mean solar time** (uniform clock time). This discrepancy arises from two effects:

1. The elliptical orbit causes varying orbital speed (eccentricity effect)
2. The tilted axis projects the sun's motion onto the equatorial plane (obliquity effect)

The equation-of-time hour angle correction [radians] is computed exactly as the difference between the mean longitude $L = M + \varpi$ and the right ascension $\alpha$ of the true Sun,

```math
\Delta \eta = L - \alpha, \qquad \alpha = \operatorname{atan2}\left(\cos\gamma \, \sin L_s,\; \cos L_s\right),
```

where $L_s = A + \varpi$ is the solar longitude. The right ascension follows from the exact projection of the ecliptic longitude onto the equatorial plane, so the obliquity contribution is valid for any tilt $\gamma$ (including $\gamma > 90°$), while the eccentricity contribution enters through the true anomaly $A$.

For small eccentricity and obliquity this reduces to the familiar perturbative form

```math
\Delta \eta \approx -2 e \sin(M) + \tan^2(\gamma/2) \sin(2M+2\varpi),
```

where the first term accounts for orbital eccentricity and the second for axial tilt. The correction can be converted to a time offset through $\Delta t = \Delta\eta T_d/(2\pi)$, where $T_d$ is the length of the solar day. For Earth it reaches roughly ±16 minutes over the year, explaining why sundials and clocks disagree.

## Hour Angle

The **hour angle** $\eta$ measures the angular distance of the sun from the local meridian (north-south line). It quantifies how far past (or before) solar noon we are at a given location:

- Local solar noon: $\eta = 0$  
- In the afternoon: $\eta > 0$
- In the morning: $\eta < 0$  

The hour angle is calculated from the time of day, with equation of time correction, and adjusted for longitude:

```math
\eta = \left( \eta_\text{uncorrected} + \Delta\eta \right) + \lambda,
```

where we have:

- The uncorrected hour angle at the prime meridian (0° longitude), $\eta_\text{uncorrected} = 2\pi t_\text{day}$ [radians]
- The fractional time of day at the prime meridian (0 at midnight, 0.5 at noon), $t_\text{day}$ [dimensionless]
- The equation of time hour angle correction, $\Delta\eta$ [radians]
- The longitude $\lambda$ [radians]

All terms are taken modulo $2\pi$ for proper angle wrapping. The factor $2\pi$ converts fractional day to angle.

## Zenith Angle

The **zenith angle** $\theta$ is the angle between the sun's rays and the vertical direction (zenith) at a location. It determines how directly sunlight strikes a surface:

- Sun directly overhead (zenith): $\theta = 0°$
- Sun at the horizon (sunrise/sunset): $\theta = 90°$
- Sun below the horizon (night): $\theta > 90°$

The cosine of the zenith angle is

```math
\cos \theta = \cos \phi \cos \delta \cos \eta + \sin \phi \sin \delta,
```

where we have:

- The latitude, $\phi$ [radians], positive northward
- The solar declination angle, $\delta$ [radians]
- The hour angle, $\eta$ [radians]

This is a fundamental equation in solar geometry. The incident solar radiation is proportional to $\cos \theta$, which is why high solar zenith angles (low sun) produce less heating than low zenith angles (overhead sun).

## Sunrise/Sunset Angle

The **sunrise/sunset hour angle** $\eta_d$ is the hour angle at which the sun crosses the horizon. It determines the length of day and night:

- Day length = $2\eta_d$ (in radians, or multiply by $T_d/(2\pi)$ for seconds)
- If $|\tan \phi \tan \delta| > 1$: polar day ($\eta_d = \pi$) or polar night ($\eta_d = 0$)

The sunrise/sunset hour angle is given by

```math
\cos \eta_d = - \tan \phi \tan \delta,
```

where this equation comes from setting $\theta = 90°$ (sun at horizon) in the zenith angle formula. The negative sign reflects that sunrise occurs at negative hour angles and sunset at positive hour angles.

## Diurnally Averaged Insolation

The **diurnally averaged** (daily-mean) insolation requires averaging $\cos \theta$ over a full 24-hour period. This is helpful for conceptual models and simpler climate models that do not resolve the full diurnal cycle.

Since insolation is proportional to $\cos \theta$, we need $\overline{\cos \theta}$, the time-averaged cosine of the zenith angle. This can be computed analytically from the sunrise/sunset hour angle:

```math
\overline{\cos \theta} = \frac{1}{\pi} \left( \eta_d \sin \phi \sin \delta + \cos \phi \cos \delta \sin \eta_d \right).
```

**Physical interpretation:**

- The $1/\pi$ prefactor combines the $1/(2\pi)$ average over the full 24-hour day with a factor of 2 from the symmetric sunrise-to-sunset integration limits ($\pm\eta_d$)
- When $\eta_d = 0$ (polar night), $\overline{\cos \theta} = 0$ (no insolation)
- When $\eta_d = \pi$ (polar day), $\overline{\cos \theta} = \sin \phi \sin \delta$ (24-hour average of the noon-and-midnight sun)
- At the equator ($\phi = 0$) during equinox ($\delta = 0$), we have $\eta_d = \pi/2$ and the formula reduces to $(1/\pi) \cos \phi \cos \delta \sin \eta_d = 1/\pi \approx 0.318$

## Azimuth Angle

The **azimuth angle** $\zeta$ specifies the horizontal direction to the sun. It is essential for tracking solar panels, understanding shading, and computing radiation on tilted surfaces. Note that this package uses a non-standard convention, measuring $\zeta$ from due East and increasing counter-clockwise (see below), rather than the more common convention of measuring clockwise from due North.

The azimuth angle is

```math
\zeta = \frac{3\pi}{2} - \operatorname{atan2}\left( \cos\delta \sin \eta,\; \cos\delta \cos \eta \sin \phi - \sin\delta \cos \phi \right),
```

where the two-argument arctangent ($\operatorname{atan2}$) resolves the correct quadrant. The numerator and denominator are written with a common factor of $\cos\delta \ge 0$ (rather than the more familiar $\sin\eta / (\cos\eta\sin\phi - \tan\delta\cos\phi)$) to avoid the $\tan\delta$ singularity at high declinations.

**Convention in this package:**

- Sun due **East**: $\zeta = 0$ (or $2\pi$)
- Sun due **North**: $\zeta = \pi/2$
- Sun due **West**: $\zeta = \pi$
- Sun due **South**: $\zeta = 3\pi/2$

The azimuth increases counter-clockwise when viewed from above. At local solar noon ($\eta = 0$), the sun is due south in the Northern Hemisphere ($\zeta = 3\pi/2$) or due north in the Southern Hemisphere.

## Planet-Star Distance

The **planet-star distance** $d$ varies throughout the year due to the elliptical orbit. This variation affects the solar flux received at the top of atmosphere through the inverse square law ($S \propto 1/d^2$).

The distance is calculated from the equation for the orbital ellipse:

```math
d = \frac{1-e^2}{1+e\cos A} d_0,
```

where we have:

- The orbital eccentricity $e$ (0 for circular, $0<e<1$ for elliptical) [unitless]
- The true anomaly $A$ [radians]
- The semi-major axis $d_0$, i.e., the mean planet-star distance

Since the insolation depends only on the *ratio* $d/d_0$ (the inverse-square law gives $S = S_0 (d_0/d)^2$ with $S_0$ the irradiance at the mean distance), `Insolation.jl` works with the dimensionless distance $d/d_0 = (1-e^2)/(1+e\cos A)$ and never needs the absolute orbit size. The returned distance is therefore expressed in units of the semi-major axis.

**For Earth:**

- Eccentricity $e \approx 0.0167$ (current value, varies over millennia)
- Perihelion (closest): $d/d_0 \approx 0.983$ (early January)
- Aphelion (farthest): $d/d_0 \approx 1.017$ (early July)
- The $\pm 1.7\%$ variation in distance causes a $\pm 3.4\%$ variation in solar flux
- The semi-major axis is $d_0 \approx 1.496 \times 10^{11}$ m $= 1$ AU; multiply $d/d_0$ by this to recover a physical distance

Earth is closest to the sun during Northern Hemisphere winter, but the obliquity effect dominates over the distance effect in determining seasons.

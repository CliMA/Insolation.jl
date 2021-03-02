# export astro_unit, planet_radius, day_length, year_anom, 
#         solar_insolation, epoch, M_epoch,
#         γ_epoch, ϖ_epoch, e_epoch
export epoch, M_epoch, γ_epoch, ϖ_epoch, e_epoch

# astro_unit() =        1.4959787e11            # 1AU [m]
# planet_radius() =     6.371e6                 # planetary radius [m]
# day_length() =        86400.0                 # length of day [s]
# year_anom() =         365.26*day_length()     # length of year [s]
# solar_insolation() =  1362.0                  # insolation [W*m-2]
epoch() =             2451545.0               # J2000 = Jan 1, 2000 at 12hr
M_epoch() =           deg2rad(357.5)          # mean anomaly at the epoch [radians]
γ_epoch() =           deg2rad(23.44)          # orbital obliquity [radians]
ϖ_epoch() =           deg2rad(282.95)         # longitude of perihelion [radians]
e_epoch() =           0.017                   # orbital eccentricity
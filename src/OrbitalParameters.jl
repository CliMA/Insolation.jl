module OrbitalParameters

export astro_unit, planet_radius, day_length, year_anom, solar_insolation

astro_unit() =        1.4959787e11            # m
planet_radius() =     6.371e6                 # m
day_length() =        86400                   # s
year_anom() =         365.26 * day_length()     # s
solar_insolation() =  1362                    # W*m-2

end
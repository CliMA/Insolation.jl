module OrbitalParameters

export astro_unit, planet_radius, day_length, year_anom, solar_insolation

"""
    Earth Parameters can be found here:
    https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
"""

astro_unit() =        1.4959787e11            # m
planet_radius() =     6.371e6                 # m
day_length() =        86400.0                 # s
year_anom() =         365.26 * day_length()   # s
solar_insolation() =  1362.0                  # W*m-2

end
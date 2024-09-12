var documenterSearchIndex = {"docs":
[{"location":"Milankovitch/#Milankovitch-Cycles","page":"Milankovitch Cycles","title":"Milankovitch Cycles","text":"","category":"section"},{"location":"Milankovitch/#Variations-in-orbital-parameters","page":"Milankovitch Cycles","title":"Variations in orbital parameters","text":"","category":"section"},{"location":"Milankovitch/","page":"Milankovitch Cycles","title":"Milankovitch Cycles","text":"using Insolation #hide\nusing Plots #hide\n\ndt = collect(-500e3:100:500e3); # years\ny = hcat(collect.(orbital_params.(dt))...);\nϖ, γ, e = y[1,:], y[2,:], y[3,:];\n\np1 = plot(dt ./ (1e3), sin.(ϖ), legend=false);\nylabel!(\"sin(ϖ)\");\np2 = plot(dt ./ (1e3), γ, legend=false);\nylabel!(\"γ\");\np3 = plot(dt ./ (1e3), e, legend=false);\nylabel!(\"e\");\nxlabel!(\"time (kyr)\")\nplot(p1, p2, p3, layout = grid(3,1), size=(600,400), dpi=250);\nsavefig(\"orbital_params.png\")","category":"page"},{"location":"Milankovitch/","page":"Milankovitch Cycles","title":"Milankovitch Cycles","text":"(Image: )","category":"page"},{"location":"Milankovitch/#Variations-in-date-of-vernal-equinox-and-perihelion-on-centennial-timescales","page":"Milankovitch Cycles","title":"Variations in date of vernal equinox and perihelion on centennial timescales","text":"","category":"section"},{"location":"Milankovitch/","page":"Milankovitch Cycles","title":"Milankovitch Cycles","text":"using Insolation #hide\nusing Plots #hide\nusing Dates #hide\nusing Roots #hide\nusing Optim #hide\n\nusing CLIMAParameters #hide\nstruct EarthParameterSet <: AbstractEarthParameterSet end #hide\nconst param_set = EarthParameterSet() #hide\n\n# Difference in NH and SH zenith angles at time x in given year\nfunction zdiff(x, year)\n    date = xtomarchdate(x,year)\n    theta_s, dist = daily_zenith_angle(date, -45., param_set, milankovitch=true)\n    theta_n, dist = daily_zenith_angle(date, 45., param_set, milankovitch=true)\n    return theta_n - theta_s\nend\n\n# x is date relative to March 1, with 1.00 representing March 1 00:00\nfunction xtomarchdate(x, year)\n    basedate = Dates.DateTime(year, 3, 1)\n    deltat = Dates.Second(round((x-1)*86400))\n    return basedate + deltat\nend\n\n# Earth-Sun distance\nfunction edist(x, year)\n    date = xtojandate(x,year)\n    _, dist = daily_zenith_angle(date, 0., param_set, milankovitch=true)\n    return dist/astro_unit()\nend\n\n# x is date relative to Jan 1, with 1.00 representing Jan 1 00:00\nfunction xtojandate(x, year)\n    basedate = Dates.DateTime(year, 1, 1)\n    deltat = Dates.Second(round((x-1)*86400))\n    date = basedate + deltat\n    return date\nend\n\nyears = 1800:2200;\ndays_eq = zeros(length(years));\ndays_per = zeros(length(years));\nfor (i,year) in enumerate(years)\n    f = (x -> zdiff(x, year))\n    days_eq[i] = find_zeros(f,-30,60)[1]\n\n    f = (x -> edist(x, year))\n    res = optimize(f,-50,50)\n    days_per[i] = Optim.minimizer(res)[1]\nend\n\nplot((years), days_eq, legend=false, dpi=250)\nxlabel!(\"Year\")\nylabel!(\"Day in March\")\ntitle!(\"Date of vernal equinox\")\nsavefig(\"equinox_dates.png\")\n\nplot((years), days_per, legend=false, dpi=250)\nxlabel!(\"Year\")\nylabel!(\"Day in Jan\")\ntitle!(\"Date of perihelion\")\nsavefig(\"perihelion_dates.png\")\n\nyears = -100e3:100:100e3 #hide\ndays_eq = zeros(length(years)) #hide\nfor (i,year) in enumerate(years) #hide\n    f = (x -> zdiff(x, year)) #hide\n    days_eq[i] = find_zeros(f,-100,100)[1] #hide\nend #hide\n\nplot((years / 1000), days_eq, legend=false, dpi=250) #hide\nxlabel!(\"kyr\") #hide\nylabel!(\"Day in March\") #hide\ntitle!(\"Date of vernal equinox\") #hide\nsavefig(\"equinox_dates_long.png\") #hide","category":"page"},{"location":"Milankovitch/","page":"Milankovitch Cycles","title":"Milankovitch Cycles","text":"(Image: ) (Image: )","category":"page"},{"location":"Milankovitch/#Variations-in-date-of-vernal-equinox-on-millenial-timescales","page":"Milankovitch Cycles","title":"Variations in date of vernal equinox on millenial timescales","text":"","category":"section"},{"location":"Milankovitch/","page":"Milankovitch Cycles","title":"Milankovitch Cycles","text":"(Image: )","category":"page"},{"location":"InsolationExamples/#Insolation-Examples","page":"Insolation Examples","title":"Insolation Examples","text":"","category":"section"},{"location":"InsolationExamples/","page":"Insolation Examples","title":"Insolation Examples","text":"These examples are based on a homework assignment from the ESE 101 class at Caltech","category":"page"},{"location":"InsolationExamples/#Insolation-in-J2000","page":"Insolation Examples","title":"Insolation in J2000","text":"","category":"section"},{"location":"InsolationExamples/","page":"Insolation Examples","title":"Insolation Examples","text":"using CLIMAParameters\nusing CLIMAParameters.Planet\nstruct EarthParameterSet <: AbstractEarthParameterSet end\nconst param_set = EarthParameterSet()\n\ninclude(\"plot_insolation.jl\")\n\nγ0 = obliq_epoch(param_set)\nϖ0 = lon_perihelion_epoch(param_set)\ne0 = eccentricity_epoch(param_set)\n\ndays, lats, F0 = calc_day_lat_insolation(365, 180, param_set)\ntitle = format(\"g = {:.2f}, w = {:.2f}, e = {:.2f}\", γ0, ϖ0, e0) #hide\nplot_day_lat_insolation(days, lats, F0, \"YlOrRd\", title, \"insol_example1.png\")","category":"page"},{"location":"InsolationExamples/","page":"Insolation Examples","title":"Insolation Examples","text":"(Image: )","category":"page"},{"location":"InsolationExamples/#Insolation-with-smaller-obliquity","page":"Insolation Examples","title":"Insolation with smaller obliquity","text":"","category":"section"},{"location":"InsolationExamples/","page":"Insolation Examples","title":"Insolation Examples","text":"using CLIMAParameters # hide\nusing CLIMAParameters.Planet # hide\nstruct EarthParameterSet <: AbstractEarthParameterSet end # hide\nconst param_set = EarthParameterSet() # hide\ninclude(\"plot_insolation.jl\") # hide\nγ0 = obliq_epoch(param_set) # hide\nϖ0 = lon_perihelion_epoch(param_set) # hide\ne0 = eccentricity_epoch(param_set) # hide\ndays, lats, F0 = calc_day_lat_insolation(365, 180, param_set) # hide\n\n# decrease γ to 20.0°\nCLIMAParameters.Planet.obliq_epoch(::EarthParameterSet) = deg2rad(20.0)\nγ1 = obliq_epoch(param_set)\ndays, lats, F2 = calc_day_lat_insolation(365, 180, param_set)\n\ntitle = format(\"g = {:.2f}, w = {:.2f}, e = {:.2f}\", γ1, ϖ0, e0) # hide\nplot_day_lat_insolation(days,lats,F2,\"YlOrRd\",  title, \"insol_example2a.png\")\ntitle = format(\"insolation diff: g' = {:.2f} - g = {:.2f}\", γ1, γ0) # hide\nplot_day_lat_insolation(days, lats, F2-F0, \"PRGn\", title, \"insol_example2b.png\")","category":"page"},{"location":"InsolationExamples/","page":"Insolation Examples","title":"Insolation Examples","text":"(Image: ) (Image: )","category":"page"},{"location":"InsolationExamples/#Insolation-with-very-large-obliquity-(like-Uranus)","page":"Insolation Examples","title":"Insolation with very large obliquity (like Uranus)","text":"","category":"section"},{"location":"InsolationExamples/","page":"Insolation Examples","title":"Insolation Examples","text":"using CLIMAParameters # hide\nusing CLIMAParameters.Planet # hide\nstruct EarthParameterSet <: AbstractEarthParameterSet end # hide\nconst param_set = EarthParameterSet() # hide\ninclude(\"plot_insolation.jl\") # hide\nγ0 = obliq_epoch(param_set) # hide\nϖ0 = lon_perihelion_epoch(param_set) # hide\ne0 = eccentricity_epoch(param_set) # hide\ndays, lats, F0 = calc_day_lat_insolation(365, 180, param_set) # hide\n\n# now change obliquity to 97.86°\nCLIMAParameters.Planet.obliq_epoch(::EarthParameterSet) = deg2rad(97.86)\nγ4 = obliq_epoch(param_set)\ndays, lats, F5 = calc_day_lat_insolation(365, 180, param_set)\n\ntitle = format(\"g = {:.2f}, w = {:.2f}, e = {:.2f}\", γ4, ϖ0, e0) # hide\nplot_day_lat_insolation(days,lats,F5,\"YlOrRd\", title, \"insol_example3a.png\")\ntitle = format(\"insolation diff: g' = {:.2f} - g = {:.2f}\", γ4, γ0) # hide\nplot_day_lat_insolation(days, lats, F5-F0, \"PRGn\", title, \"insol_example3b.png\")\n","category":"page"},{"location":"InsolationExamples/","page":"Insolation Examples","title":"Insolation Examples","text":"(Image: ) (Image: )","category":"page"},{"location":"library/#Application-Programming-Interface-(APIs)","page":"APIs","title":"Application Programming Interface (APIs)","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Documenting the public user interface","category":"page"},{"location":"library/#Orbital-Parameters","page":"APIs","title":"Orbital Parameters","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [Insolation]\nPrivate = false\nPages   = [\"Insolation.jl\"]","category":"page"},{"location":"library/#Insolation.orbital_params-Tuple{FT} where FT<:Real","page":"APIs","title":"Insolation.orbital_params","text":"orbital_params(dt::FT) where {FT <: Real}\n\nThis function returns the orbital parameters (ϖ, γ, e) at dt (years) since J2000 epoch. Data are read from file and interpolation function are created in init() method. The functions are stored as global variables that are used inside of Insolation.jl. The parameters vary due to Milankovitch cycles.  Data from this paper are in the \"src/data/INSOL.LA2004.BTL.csv\" file.\n\n\n\n\n\n","category":"method"},{"location":"library/#Zenith-Angle","page":"APIs","title":"Zenith Angle","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [Insolation]\nPrivate = false\nPages   = [\"ZenithAngleCalc.jl\"]","category":"page"},{"location":"library/#Insolation.daily_zenith_angle-Union{Tuple{FT}, Tuple{Dates.DateTime,FT,CLIMAParameters.AbstractParameterSet}} where FT<:Real","page":"APIs","title":"Insolation.daily_zenith_angle","text":"daily_zenith_angle(date::DateTime,\n                   latitude::FT,\n                   param_set::APS;\n                   eot_correction::Bool=true,\n                   milankovitch::Bool=true) where {FT <: Real}\n\nReturns the daily averaged zenith angle and earth-sun distance at a particular latitude given the date and orbital parameters obliquity, longitude of perihelion, and eccentricity param_set is an AbstractParameterSet from CLIMAParameters.jl.\n\neot_correction is an optional Boolean keyword argument that defaults to true when set to true the equation of time correction is turned on. This switch functionality is implemented for easy comparisons with reanalyses.\n\nmilankovitch is an optional Boolean keyword argument that defaults to true when set to true the orbital parameters are calculated for the given DateTime, when set to false the orbital parameters at the J2000 epoch from CLIMAParameters are used.\n\n\n\n\n\n","category":"method"},{"location":"library/#Insolation.instantaneous_zenith_angle-Union{Tuple{FT}, Tuple{Dates.DateTime,FT,FT,CLIMAParameters.AbstractParameterSet}} where FT<:Real","page":"APIs","title":"Insolation.instantaneous_zenith_angle","text":"instantaneous_zenith_angle(date::DateTime,\n                           longitude::FT,\n                           latitude::FT,\n                           param_set::APS;\n                           eot_correction::Bool=true,\n                           milankovitch::Bool=true) where {FT <: Real}\n\nReturns the zenith angle and earth-sun distance at a particular longitude and latitude on the given date (and time UTC) given orbital parameters: obliquity, longitude of perihelion, and eccentricity param_set is an AbstractParameterSet from CLIMAParameters.jl.\n\neot_correction is an optional Boolean keyword argument that defaults to true when set to true the equation of time correction is turned on. This switch functionality is implemented for easy comparisons with reanalyses.\n\nmilankovitch is an optional Boolean keyword argument that defaults to true when set to true the orbital parameters are calculated for the given DateTime when set to false the orbital parameters at the J2000 epoch from CLIMAParameters are used.\n\n\n\n\n\n","category":"method"},{"location":"library/#Insolation","page":"APIs","title":"Insolation","text":"","category":"section"},{"location":"library/","page":"APIs","title":"APIs","text":"Modules = [Insolation]\nPrivate = false\nPages   = [\"InsolationCalc.jl\"]","category":"page"},{"location":"library/#Insolation.insolation-Union{Tuple{FT}, Tuple{FT,FT,CLIMAParameters.AbstractParameterSet}} where FT<:Real","page":"APIs","title":"Insolation.insolation","text":"insolation(θ::FT, d::FT, param_set::APS) where {FT <: Real}\n\nReturns the insolation given the zenith angle and earth-sun distance param_set is an AbstractParameterSet from CLIMAParameters.jl.\n\n\n\n\n\n","category":"method"},{"location":"library/#Insolation.solar_flux_and_cos_sza-Union{Tuple{FT}, Tuple{Dates.DateTime,FT,FT,CLIMAParameters.AbstractParameterSet}} where FT<:Real","page":"APIs","title":"Insolation.solar_flux_and_cos_sza","text":"solar_flux_and_cos_sza(date::DateTime,\n                  longitude::FT,\n                  latitude::FT,\n                  param_set::APS) where {FT <: Real}\n\nReturns the top-of-atmosphere (TOA) solar flux, i.e.  the total solar irradiance (TSI) weighted by the earth-sun distance and cos(solar zenith angle) for input to RRTMGP.jl param_set is an AbstractParameterSet from CLIMAParameters.jl.\n\n\n\n\n\n","category":"method"},{"location":"ZenithAngleEquations/#Zenith-Angle-Equations","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"","category":"section"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"These equations are from Tapio Schneiders's textbook draft chapter 3.","category":"page"},{"location":"ZenithAngleEquations/#Mean-Anomaly-(3.6)","page":"Zenith Angle Equations","title":"Mean Anomaly (3.6)","text":"","category":"section"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"The mean anomaly M at current time t is,","category":"page"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"M = frac2pi (t - t_0)Y_a + M_0","category":"page"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"where t_0 is the time at the epoch (J2000), defined as January 1, 2000 at 12hr UTC,  M_0 is the mean anomaly at the epoch, and Y_a is the length of the anomalistic year.","category":"page"},{"location":"ZenithAngleEquations/#True-Anomaly-(3.8)","page":"Zenith Angle Equations","title":"True Anomaly (3.8)","text":"","category":"section"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"The true anomaly A is given by,","category":"page"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"A = M + left( 2e - frac14e^3 right) sin(M) + frac54 e^2 sin(2M) + frac1312 e^3 sin(3M)","category":"page"},{"location":"ZenithAngleEquations/#True-Longitude-(3.9)","page":"Zenith Angle Equations","title":"True Longitude (3.9)","text":"","category":"section"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"The true longitude is the sum of the true anomaly and the longitude of perihelion (varpi).","category":"page"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"L = A + varpi","category":"page"},{"location":"ZenithAngleEquations/#Declination-(3.16)","page":"Zenith Angle Equations","title":"Declination (3.16)","text":"","category":"section"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"The sin of the declination angle is,","category":"page"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"sin delta = sin gamma sin L","category":"page"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"where gamma is the orbital obliquity.","category":"page"},{"location":"ZenithAngleEquations/#Equation-of-Time","page":"Zenith Angle Equations","title":"Equation of Time","text":"","category":"section"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"The equation of time corrects for the difference between apparent solar time (local solar noon) and mean solar time due to the obliquity of the planet and eccentricity of the orbit.","category":"page"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"Delta t = -2 e sin(M) + tan^2(gamma2) sin(2M+2varpi)","category":"page"},{"location":"ZenithAngleEquations/#Hour-Angle-(3.17)","page":"Zenith Angle Equations","title":"Hour Angle (3.17)","text":"","category":"section"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"The hour angle eta defines how the sun's position changes during a day due to the planet's rotation.","category":"page"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"eta = frac2pi (t+Delta t)T_d","category":"page"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"where t is the current time referenced to the epoch t_0 and T_d is the length of a day.","category":"page"},{"location":"ZenithAngleEquations/#Zenith-Angle-(3.18)","page":"Zenith Angle Equations","title":"Zenith Angle (3.18)","text":"","category":"section"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"The zenith angle at any time is given by,","category":"page"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"cos theta = cos phi cos delta cos eta + sin phi sin delta","category":"page"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"where phi is the latitude.","category":"page"},{"location":"ZenithAngleEquations/#Sunrise/Sunset-Angle-(3.19)","page":"Zenith Angle Equations","title":"Sunrise/Sunset Angle (3.19)","text":"","category":"section"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"The sunrise/sunset angle is the hour angle eta_d at which the sun rises or sets,","category":"page"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"cos eta_d = - tan phi tan delta","category":"page"},{"location":"ZenithAngleEquations/#Daily-averaged-Zenith-Angle-(3.20)","page":"Zenith Angle Equations","title":"Daily-averaged Zenith Angle (3.20)","text":"","category":"section"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"The daily-averaged zenith angle can be defined then in terms of the sunrise/sunset angle as,","category":"page"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"overlinecos theta = frac1pi left( eta_d sin phi sin delta + cos phi cos delta cos eta_d right)","category":"page"},{"location":"ZenithAngleEquations/#Azimuth-Angle","page":"Zenith Angle Equations","title":"Azimuth Angle","text":"","category":"section"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"The azimuth angle is,","category":"page"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"zeta = frac3pi2 - arctan left( fracsin etacos eta sin phi - tan delta cos phi right)","category":"page"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"The azimuth is defined as 0 to the East and increasing counter-clockwise, such that at local solar noon when eta=0, then zeta = frac3pi2.","category":"page"},{"location":"ZenithAngleEquations/#Planet-Sun-Distance-(3.1)","page":"Zenith Angle Equations","title":"Planet-Sun Distance (3.1)","text":"","category":"section"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"The distance between the planet and the sun, referenced to the mean distance, as defined by the eccentricity of the orbit.","category":"page"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"d = frac1-e^21+ecos A d_0","category":"page"},{"location":"ZenithAngleEquations/","page":"Zenith Angle Equations","title":"Zenith Angle Equations","text":"For the Earth, d_0 = 1 AU.","category":"page"},{"location":"#Insolation.jl","page":"Home","title":"Insolation.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Insolation","category":"page"},{"location":"","page":"Home","title":"Home","text":"Insolation.jl is a library that calculates the zenith angle and insolation  at a given point on Earth (lat, lon) for a given date/time and orbital configuration  (obliquity, eccentricity, and longitude of perihelion). The library is split between two files, ZenithAngleCalc.jl which calculates the zenith angle, azimuth angle, and Earth-sun distance,  and InsolationCalc.jl which calculates the insolation given a zenith angle.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The zenith angle and insolation can both be calculated either as instantaneous  values or as daily averaged values. The functions in ZenithAngleCalc.jl are  overwritten to accept a variety of inputs, either calculating the orbital parameters  given a DateTime object, or prescribing the orbital parameters as input.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The equations used for calculating orbital parameters are from Tapio Schneider's textbook draft.  See Zenith Angle Equations for more details.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package has been designed to be flexible for non-Earth settings as well. This specific functionality has not yet been implemented [as of 13 Sept 2021].","category":"page"},{"location":"#Authors","page":"Home","title":"Authors","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Insolation.jl is being developed by the Climate Modeling Alliance. Specifically it has been developed to be used by RRTMGP.jl  and the ClimateMachine.jl.","category":"page"}]
}
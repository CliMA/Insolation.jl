var documenterSearchIndex = {"docs":
[{"location":"InsolationExamples/#Insolation-Examples-1","page":"Insolation Examples","title":"Insolation Examples","text":"","category":"section"},{"location":"InsolationExamples/#Example-1-1","page":"Insolation Examples","title":"Example 1","text":"","category":"section"},{"location":"InsolationExamples/#","page":"Insolation Examples","title":"Insolation Examples","text":"include(\"plot_insolation.jl\")\n\n# 2018 orbital constants\nγ0 = 23.44\nϖ0 = 282.95\ne0 = 0.017\n\ndays, lats, F0 = calc_day_lat_insolation(365, 180, γ0, ϖ0, e0)\ntitle = \"obliq=\" * \"$(γ0)\" * \", perihelion=\" * \"$(ϖ0)\" * \", ecc=\" * \"$(e0)\"\nplot_day_lat_insolation(days, lats, F0, \"YlOrRd\", title, \"insol_example1.png\")","category":"page"},{"location":"InsolationExamples/#","page":"Insolation Examples","title":"Insolation Examples","text":"(Image: )","category":"page"},{"location":"InsolationExamples/#Example-2-1","page":"Insolation Examples","title":"Example 2","text":"","category":"section"},{"location":"InsolationExamples/#","page":"Insolation Examples","title":"Insolation Examples","text":"include(\"plot_insolation.jl\") # hide\n\n# turn longitude of perihelion by 180°\nγ0 = 23.44\nϖ0 = 282.95 # hide\nϖ1 = ϖ0 + 180.0\ne0 = 0.017\n\ndays, lats, F0 = calc_day_lat_insolation(365, 180, γ0, ϖ0, e0) # hide\ndays, lats, F1 = calc_day_lat_insolation(365, 180, γ0, ϖ1, e0)\n\ntitle = \"obliq=\" * \"$(γ0)\" * \", perihelion=\" * \"$(ϖ1)\" * \", ecc=\" * \"$(e0)\"\nplot_day_lat_insolation(days, lats, F1, \"YlOrRd\",  title, \"insol_example2a.png\")\n\ntitle = \"insolation diff: perihelion0=\" * \"$(ϖ0)\" * \", perihelion1=\" * \"$(ϖ1)\"\nplot_day_lat_insolation(days,lats,F1-F0,\"PRGn\", title, \"insol_example2b.png\")","category":"page"},{"location":"InsolationExamples/#","page":"Insolation Examples","title":"Insolation Examples","text":"(Image: ) (Image: )","category":"page"},{"location":"InsolationExamples/#Example-3-1","page":"Insolation Examples","title":"Example 3","text":"","category":"section"},{"location":"InsolationExamples/#","page":"Insolation Examples","title":"Insolation Examples","text":"include(\"plot_insolation.jl\") # hide\n\n# perihelion back to normal. decrease γ to 22.0°\nγ0 = 23.44 # hide\nγ1 = 22.0\nϖ0 = 282.95\ne0 = 0.017\n\ndays, lats, F0 = calc_day_lat_insolation(365, 180, γ0, ϖ0, e0) # hide\ndays, lats, F2 = calc_day_lat_insolation(365, 180, γ1, ϖ0, e0)\n\ntitle = \"obliq=\" * \"$(γ1)\" * \", perihelion=\" * \"$(ϖ0)\" * \", ecc=\" * \"$(e0)\"\nplot_day_lat_insolation(days,lats,F2,\"YlOrRd\",  title, \"insol_example3a.png\")\n\ntitle = \"insolation diff: obliq0=\" * \"$(γ0)\" * \", obliq1=\" * \"$(γ1)\"\nplot_day_lat_insolation(days,lats,F2-F0,\"PRGn\", title, \"insol_example3b.png\")","category":"page"},{"location":"InsolationExamples/#","page":"Insolation Examples","title":"Insolation Examples","text":"(Image: ) (Image: )","category":"page"},{"location":"InsolationExamples/#Example-4-1","page":"Insolation Examples","title":"Example 4","text":"","category":"section"},{"location":"InsolationExamples/#","page":"Insolation Examples","title":"Insolation Examples","text":"include(\"plot_insolation.jl\") # hide\n\n# decrease γ further to 18.0°\nγ0 = 23.44 # hide\nγ2 = 18.0\nϖ0 = 282.95\ne0 = 0.017\n\ndays, lats, F0 = calc_day_lat_insolation(365, 180, γ0, ϖ0, e0) # hide\ndays, lats, F3 = calc_day_lat_insolation(365, 180, γ2, ϖ0, e0)\n\ntitle = \"obliq=\" * \"$(γ2)\" * \", perihelion=\" * \"$(ϖ0)\" * \", ecc=\" * \"$(e0)\"\nplot_day_lat_insolation(days,lats,F3,\"YlOrRd\", title, \"insol_example4a.png\")\n\ntitle = \"insolation diff: obliq0=\" * \"$(γ0)\" * \", obliq1=\" * \"$(γ2)\"\nplot_day_lat_insolation(days,lats,F3-F0,\"PRGn\", title, \"insol_example4b.png\")","category":"page"},{"location":"InsolationExamples/#","page":"Insolation Examples","title":"Insolation Examples","text":"(Image: ) (Image: )","category":"page"},{"location":"InsolationExamples/#Example-5-1","page":"Insolation Examples","title":"Example 5","text":"","category":"section"},{"location":"InsolationExamples/#","page":"Insolation Examples","title":"Insolation Examples","text":"include(\"plot_insolation.jl\") # hide\n\n# now change obliquity to 60.0°\nγ0 = 23.44 # hide\nγ3 = 60.0\nϖ0 = 282.95\ne0 = 0.017\n\ndays, lats, F0 = calc_day_lat_insolation(365, 180, γ0, ϖ0, e0) # hide\ndays, lats, F4 = calc_day_lat_insolation(365, 180, γ3, ϖ0, e0)\n\ntitle = \"obliq=\" * \"$(γ3)\" * \", perihelion=\" * \"$(ϖ0)\" * \", ecc=\" * \"$(e0)\"\nplot_day_lat_insolation(days,lats,F4,\"YlOrRd\", title, \"insol_example5a.png\")\n\ntitle = \"insolation diff: obliq0=\" * \"$(γ0)\" * \", obliq1=\" * \"$(γ3)\"\nplot_day_lat_insolation(days,lats,F4-F0,\"PRGn\", title, \"insol_example5b.png\")","category":"page"},{"location":"InsolationExamples/#","page":"Insolation Examples","title":"Insolation Examples","text":"(Image: ) (Image: )","category":"page"},{"location":"InsolationExamples/#Example-6-1","page":"Insolation Examples","title":"Example 6","text":"","category":"section"},{"location":"InsolationExamples/#","page":"Insolation Examples","title":"Insolation Examples","text":"include(\"plot_insolation.jl\") # hide\n\n# now change obliquity to 97.86°\nγ0 = 23.44 # hide\nγ4 = 97.86\nϖ0 = 282.95\ne0 = 0.017\n\ndays, lats, F0 = calc_day_lat_insolation(365, 180, γ0, ϖ0, e0) # hide\ndays, lats, F5 = calc_day_lat_insolation(365, 180, γ4, ϖ0, e0)\n\ntitle = \"obliq=\" * \"$(γ4)\" * \", perihelion=\" * \"$(ϖ0)\" * \", ecc=\" * \"$(e0)\"\nplot_day_lat_insolation(days,lats,F5,\"YlOrRd\", title, \"insol_example6a.png\")\n\ntitle = \"insolation diff: obliq0=\" * \"$(γ0)\" * \", obliq1=\" * \"$(γ4)\"\nplot_day_lat_insolation(days,lats,F5-F0,\"PRGn\", title, \"insol_example6b.png\")\n","category":"page"},{"location":"InsolationExamples/#","page":"Insolation Examples","title":"Insolation Examples","text":"(Image: ) (Image: )","category":"page"},{"location":"library/#Library-1","page":"Library","title":"Library","text":"","category":"section"},{"location":"library/#","page":"Library","title":"Library","text":"Documenting the public user interface","category":"page"},{"location":"library/#Solar-Zenith-Angle-1","page":"Library","title":"Solar Zenith Angle","text":"","category":"section"},{"location":"library/#","page":"Library","title":"Library","text":"Modules = [Insolation.SolarZenithAngle]\nPrivate = false\nPages   = [\"SolarZenithAngle.jl\"]","category":"page"},{"location":"library/#Insolation.SolarZenithAngle.daily_zenith_angle-Union{Tuple{FT}, Tuple{Dates.DateTime,FT}} where FT<:Real","page":"Library","title":"Insolation.SolarZenithAngle.daily_zenith_angle","text":"daily_zenith_angle(date::DateTime,\n                   latitude::FT) where {FT <: Real}\n\nreturns the daily averaged zenith angle and earth-sun distance at a particular latitude on the given date\n\n\n\n\n\n","category":"method"},{"location":"library/#Insolation.SolarZenithAngle.daily_zenith_angle-Union{Tuple{I}, Tuple{FT}, Tuple{I,FT,FT,FT,FT}} where I<:Int64 where FT<:Real","page":"Library","title":"Insolation.SolarZenithAngle.daily_zenith_angle","text":"daily_zenith_angle(days_since_equinox::I,\n                   obliquity::FT,\n                   perihelion::FT,\n                   eccentricity::FT,\n                   latitude::FT) where {FT <: Real, I <: Int}\n\nreturns the daily averaged zenith angle and earth-sun distance at a particular latitude given the days since vernal equinox (defined as March 21), orbital obliquity, longitude of perihelion, and eccentricity\n\n\n\n\n\n","category":"method"},{"location":"library/#Insolation.SolarZenithAngle.instantaneous_zenith_angle-Union{Tuple{FT}, Tuple{Dates.DateTime,FT,FT,FT,FT,FT,FT}} where FT<:Real","page":"Library","title":"Insolation.SolarZenithAngle.instantaneous_zenith_angle","text":"instantaneous_zenith_angle(date::DateTime,\n                           timezone::FT,\n                           longitude::FT,\n                           latitude::FT,\n                           obliquity::FT,\n                           perihelion::FT,\n                           eccentricity::FT) where {FT <: Real}\n\nreturns the zenith angle and earth-sun distance at a particular longitude and latitude on the given date given orbital parameters: obliquity, longitude of perihelion, and eccentricity\n\n\n\n\n\n","category":"method"},{"location":"library/#Insolation.SolarZenithAngle.instantaneous_zenith_angle-Union{Tuple{FT}, Tuple{Dates.DateTime,FT,FT,FT}} where FT<:Real","page":"Library","title":"Insolation.SolarZenithAngle.instantaneous_zenith_angle","text":"instantaneous_zenith_angle(date::DateTime,\n                           timezone::FT,\n                           longitude::FT,\n                           latitude::FT) where {FT <: Real}\n\nreturns the zenith angle and earth-sun distance at a particular longitude and latitude on the given date\n\nadd citations: \n\nhttps://www.esrl.noaa.gov/gmd/grad/solcalc/\nhttps://www.cfa.harvard.edu/~jzhao/times.html\nhttps://github.com/thabbott/zenithangle/blob/master/solar.js\nhttps://github.com/claresinger/3d-cloud-rad/blob/master/sza-scripts/sza_utils.py\nhttp://farside.ph.utexas.edu/Books/Syntaxis/Almagest/node36.html\nhttps://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html\n\"Astronomical Algorithms\" by Jean Meeus\n\n\n\n\n\n","category":"method"},{"location":"library/#Solar-Insolation-1","page":"Library","title":"Solar Insolation","text":"","category":"section"},{"location":"library/#","page":"Library","title":"Library","text":"Modules = [Insolation.SolarInsolation]\nPrivate = false\nPages   = [\"SolarInsolation.jl\"]","category":"page"},{"location":"library/#Insolation.SolarInsolation.daily_insolation-Union{Tuple{FT}, Tuple{Dates.DateTime,FT}} where FT<:Real","page":"Library","title":"Insolation.SolarInsolation.daily_insolation","text":"daily_insolation(date::DateTime,\n                 latitude::FT) where {FT <: Real}\n\nreturns the daily averaged insolation at a particular latitude on the given date\n\n\n\n\n\n","category":"method"},{"location":"library/#Insolation.SolarInsolation.daily_insolation-Union{Tuple{I}, Tuple{FT}, Tuple{I,FT,FT,FT,FT}} where I<:Int64 where FT<:Real","page":"Library","title":"Insolation.SolarInsolation.daily_insolation","text":"daily_insolation(days_since_equinox::I,\n                 obliquity::FT,\n                 perihelion::FT,\n                 eccentricity::FT,\n                 latitude::FT) where {FT <: Real, I <: Int}\n\nreturns the daily averaged insolation at a particular latitude given the days since vernal equinox (defined as March 21), orbital obliquity, longitude of perihelion, and eccentricity\n\n\n\n\n\n","category":"method"},{"location":"library/#Insolation.SolarInsolation.insolation-Union{Tuple{FT}, Tuple{FT,FT}} where FT<:Real","page":"Library","title":"Insolation.SolarInsolation.insolation","text":"insolation(θ::FT,\n                 d::FT) where {FT <: Real}\n\nreturns the insolation given the zenith angle and earth-sun distance\n\n\n\n\n\n","category":"method"},{"location":"library/#Insolation.SolarInsolation.instantaneous_insolation-Union{Tuple{FT}, Tuple{Dates.DateTime,FT,FT,FT}} where FT<:Real","page":"Library","title":"Insolation.SolarInsolation.instantaneous_insolation","text":"instantaneous_insolation(date::DateTime,\n                         timezone::FT,\n                         longitude::FT,\n                         latitude::FT) where {FT <: Real}\n\nreturns the daily averaged insolation at a particular latitude on the given date\n\n\n\n\n\n","category":"method"},{"location":"ZenithAngleExamples/#Solar-Zenith-Angle-Examples-1","page":"Zenith Angle Examples","title":"Solar Zenith Angle Examples","text":"","category":"section"},{"location":"ZenithAngleExamples/#Example-1-1","page":"Zenith Angle Examples","title":"Example 1","text":"","category":"section"},{"location":"ZenithAngleExamples/#","page":"Zenith Angle Examples","title":"Zenith Angle Examples","text":"include(\"plot_zenith.jl\")\n\n# 2018 orbital constants\nγ0 = 23.44\nϖ0 = 282.95\ne0 = 0.017\n\ndays, lats, S = calc_day_lat_sza(365, 180, γ0, ϖ0, e0)\ntitle = \"obliq=\" * \"$(γ0)\" * \", perihelion=\" * \"$(ϖ0)\" * \", ecc=\" * \"$(e0)\"\nplot_day_lat_sza(days, lats, S, \"YlOrRd\", title, \"sza_example1.png\")","category":"page"},{"location":"ZenithAngleExamples/#","page":"Zenith Angle Examples","title":"Zenith Angle Examples","text":"(Image: )","category":"page"},{"location":"ZenithAngleExamples/#Example-2-1","page":"Zenith Angle Examples","title":"Example 2","text":"","category":"section"},{"location":"ZenithAngleExamples/#","page":"Zenith Angle Examples","title":"Zenith Angle Examples","text":"include(\"plot_zenith.jl\") # hide\n\n# turn longitude of perihelion by 180°\nγ0 = 23.44\nϖ0 = 282.95 # hide\nϖ1 = ϖ0 + 180.0\ne0 = 0.017\n\ndays, lats, S0 = calc_day_lat_sza(365, 180, γ0, ϖ0, e0) # hide\ndays, lats, S1 = calc_day_lat_sza(365, 180, γ0, ϖ1, e0)\n\ntitle = \"obliq=\" * \"$(γ0)\" * \", perihelion=\" * \"$(ϖ1)\" * \", ecc=\" * \"$(e0)\"\nplot_day_lat_sza(days, lats, S1, \"YlOrRd\",  title, \"sza_example2a.png\")\n\ntitle = \"SZA diff: perihelion0=\" * \"$(ϖ0)\" * \", perihelion1=\" * \"$(ϖ1)\"\nplot_day_lat_sza(days,lats,S1-S0,\"PRGn\", title, \"sza_example2b.png\")","category":"page"},{"location":"ZenithAngleExamples/#","page":"Zenith Angle Examples","title":"Zenith Angle Examples","text":"(Image: ) (Image: )","category":"page"},{"location":"ZenithAngleExamples/#Example-3-1","page":"Zenith Angle Examples","title":"Example 3","text":"","category":"section"},{"location":"ZenithAngleExamples/#","page":"Zenith Angle Examples","title":"Zenith Angle Examples","text":"include(\"plot_zenith.jl\") # hide\n\n# perihelion back to normal. decrease γ to 22.0°\nγ0 = 23.44 # hide\nγ1 = 22.0\nϖ0 = 282.95\ne0 = 0.017\n\ndays, lats, S0 = calc_day_lat_sza(365, 180, γ0, ϖ0, e0) # hide\ndays, lats, S2 = calc_day_lat_sza(365, 180, γ1, ϖ0, e0)\n\ntitle = \"obliq=\" * \"$(γ1)\" * \", perihelion=\" * \"$(ϖ0)\" * \", ecc=\" * \"$(e0)\"\nplot_day_lat_sza(days,lats,S2,\"YlOrRd\",  title, \"sza_example3a.png\")\n\ntitle = \"SZA diff: obliq0=\" * \"$(γ0)\" * \", obliq1=\" * \"$(γ1)\"\nplot_day_lat_sza(days,lats,S2-S0,\"PRGn\", title, \"sza_example3b.png\")","category":"page"},{"location":"ZenithAngleExamples/#","page":"Zenith Angle Examples","title":"Zenith Angle Examples","text":"(Image: ) (Image: )","category":"page"},{"location":"ZenithAngleExamples/#Example-4-1","page":"Zenith Angle Examples","title":"Example 4","text":"","category":"section"},{"location":"ZenithAngleExamples/#","page":"Zenith Angle Examples","title":"Zenith Angle Examples","text":"include(\"plot_zenith.jl\") # hide\n\n# decrease γ further to 18.0°\nγ0 = 23.44 # hide\nγ2 = 18.0\nϖ0 = 282.95\ne0 = 0.017\n\ndays, lats, S0 = calc_day_lat_sza(365, 180, γ0, ϖ0, e0) # hide\ndays, lats, S3 = calc_day_lat_sza(365, 180, γ2, ϖ0, e0)\n\ntitle = \"obliq=\" * \"$(γ2)\" * \", perihelion=\" * \"$(ϖ0)\" * \", ecc=\" * \"$(e0)\"\nplot_day_lat_sza(days,lats,S3,\"YlOrRd\", title, \"sza_example4a.png\")\n\ntitle = \"SZA diff: obliq0=\" * \"$(γ0)\" * \", obliq1=\" * \"$(γ2)\"\nplot_day_lat_sza(days,lats,S3-S0,\"PRGn\", title, \"sza_example4b.png\")","category":"page"},{"location":"ZenithAngleExamples/#","page":"Zenith Angle Examples","title":"Zenith Angle Examples","text":"(Image: ) (Image: )","category":"page"},{"location":"ZenithAngleExamples/#Example-5-1","page":"Zenith Angle Examples","title":"Example 5","text":"","category":"section"},{"location":"ZenithAngleExamples/#","page":"Zenith Angle Examples","title":"Zenith Angle Examples","text":"include(\"plot_zenith.jl\") # hide\n\n# now change obliquity to 60.0°\nγ0 = 23.44 # hide\nγ3 = 60.0\nϖ0 = 282.95\ne0 = 0.017\n\ndays, lats, S0 = calc_day_lat_sza(365, 180, γ0, ϖ0, e0) # hide\ndays, lats, S4 = calc_day_lat_sza(365, 180, γ3, ϖ0, e0)\n\ntitle = \"obliq=\" * \"$(γ3)\" * \", perihelion=\" * \"$(ϖ0)\" * \", ecc=\" * \"$(e0)\"\nplot_day_lat_sza(days,lats,S4,\"YlOrRd\", title, \"sza_example5a.png\")\n\ntitle = \"SZA diff: obliq0=\" * \"$(γ0)\" * \", obliq1=\" * \"$(γ3)\"\nplot_day_lat_sza(days,lats,S4-S0,\"PRGn\", title, \"sza_example5b.png\")","category":"page"},{"location":"ZenithAngleExamples/#","page":"Zenith Angle Examples","title":"Zenith Angle Examples","text":"(Image: ) (Image: )","category":"page"},{"location":"ZenithAngleExamples/#Example-6-1","page":"Zenith Angle Examples","title":"Example 6","text":"","category":"section"},{"location":"ZenithAngleExamples/#","page":"Zenith Angle Examples","title":"Zenith Angle Examples","text":"include(\"plot_zenith.jl\") # hide\n\n# now change obliquity to 97.86°\nγ0 = 23.44 # hide\nγ4 = 97.86\nϖ0 = 282.95\ne0 = 0.017\n\ndays, lats, S0 = calc_day_lat_sza(365, 180, γ0, ϖ0, e0) # hide\ndays, lats, S5 = calc_day_lat_sza(365, 180, γ4, ϖ0, e0)\n\ntitle = \"obliq=\" * \"$(γ4)\" * \", perihelion=\" * \"$(ϖ0)\" * \", ecc=\" * \"$(e0)\"\nplot_day_lat_sza(days,lats,S5,\"YlOrRd\", title, \"sza_example6a.png\")\n\ntitle = \"SZA diff: obliq0=\" * \"$(γ0)\" * \", obliq1=\" * \"$(γ4)\"\nplot_day_lat_sza(days,lats,S5-S0,\"PRGn\", title, \"sza_example6b.png\")\n","category":"page"},{"location":"ZenithAngleExamples/#","page":"Zenith Angle Examples","title":"Zenith Angle Examples","text":"(Image: ) (Image: )","category":"page"},{"location":"#Insolation.jl-1","page":"Home","title":"Insolation.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"CurrentModule = Insolation","category":"page"}]
}
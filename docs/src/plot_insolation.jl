using Insolation.ZenithAngleCalc
using Insolation.InsolationCalc
using Insolation.OrbitalParameters
using Plots
using Dates
using Formatting

"""
    calc_day_lat_insolation(n_days::I,
                            n_lats::I,
                            γ::FT,
                            ϖ::FT,
                            e::FT) where {FT<:Real,I<:Int}
"""
function calc_day_lat_insolation(n_days::I,
                                n_lats::I,
                                γ::FT,
                                ϖ::FT,
                                e::FT) where {FT<:Real,I<:Int}
  d_arr = Array{I}(round.(collect(range(0, stop = 365, length = n_days))))
  l_arr = Array{I}(round.(collect(range(-90, stop = 90, length = n_lats))))
  l_arr = collect(range(-90, stop = 90, length = n_lats))
  F_arr = zeros(FT, n_days, n_lats)
  # loop over days
  for (i, d) in enumerate(d_arr)
    for (j, lat) in enumerate(l_arr)
        date = Dates.DateTime(2020,1,1) + Dates.Day(d)
        θ, dist = daily_zenith_angle(date, lat, γ, ϖ, e)
        F_arr[i, j] = insolation(θ, dist)
    end
  end
  return d_arr, l_arr, F_arr
end

"""
    plot_day_lat_insolation(n_days::I,
                            n_lats::I,
                            F_arr::Array{FT},
                            scmap,
                            stitle,
                            file_name) where {FT<:AbstractFloat,I<:Int}
"""
function plot_day_lat_insolation(d_arr::Array{I},
                                 l_arr::Array{FT},
                                 F_arr::Array{FT},
                                 scmap,
                                 stitle,
                                 file_name) where {FT<:AbstractFloat,I<:Int}
  if scmap == "YlOrRd"
    cmap = :YlOrRd
    vmin, vmax = 0, ceil(max(F_arr...)/100)*100
  elseif scmap == "PRGn"
    cmap = :PRGn
    vmin, vmax = ceil(max(abs.(F_arr)...)/10)*-10, ceil(max(abs.(F_arr)...)/10)*10
  end
  
  p = contourf(d_arr, l_arr, F_arr', c=cmap, clims=(vmin,vmax), title=stitle, 
    xlabel="Days since Jan 1 2020", ylabel="Latitude", colorbar_title="ToA Insolation [W/m2]")

  savefig(file_name)
end

# function main()
#   γ0 = obliquity()
#   ϖ0 = perihelion()
#   e0 = eccentricity()
#   days, lats, F0 = calc_day_lat_insolation(365, 180, γ0, ϖ0, e0)
#   title = format("g = {:.2f}, w = {:.2f}, e = {:.2f}", γ0, ϖ0, e0)
#   plot_day_lat_insolation(days, lats, F0, "YlOrRd", title, "insol_example1.png")

#   ϖ1 = perihelion() + π
#   days, lats, F1 = calc_day_lat_insolation(365, 180, γ0, ϖ1, e0)
#   title = format("g = {:.2f}, w = {:.2f}, e = {:.2f}", γ0, ϖ1, e0)
#   plot_day_lat_insolation(days, lats, F1, "YlOrRd",  title, "insol_example2a.png")

#   title = format("insolation diff: w0 = {:.2f}, w1 = {:.2f}", ϖ0, ϖ1)
#   plot_day_lat_insolation(days,lats,F1-F0,"PRGn", title, "insol_example2b.png")
# end

# main()
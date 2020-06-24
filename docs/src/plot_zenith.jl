using Insolation.SolarZenithAngle
using Plots

"""
    calc_day_lat_sza(n_days::I,
                            n_lats::I,
                            γ::FT,
                            ϖ::FT,
                            e::FT) where {FT<:AbstractFloat,I<:Int}
"""
function calc_day_lat_sza(n_days::I,
                                n_lats::I,
                                γ::FT,
                                ϖ::FT,
                                e::FT) where {FT<:AbstractFloat,I<:Int}
  equinox_day = 76
  d_arr = Array{I}(round.(collect(range(0, stop = 365, length = n_days))))
  l_arr = Array{I}(round.(collect(range(-90, stop = 90, length = n_lats))))
  l_arr = collect(range(-90, stop = 90, length = n_lats))
  sza_arr = zeros(FT, n_days, n_lats)
  # loop over days
  for (i, day) in enumerate(d_arr)
    for (j, lat) in enumerate(l_arr)
      sza_rad = daily_zenith_angle(day-equinox_day, γ, ϖ, e, lat)[1]
      sza_arr[i, j] = rad2deg(sza_rad)
    end
  end
  return d_arr, l_arr, sza_arr
end

"""
    plot_day_lat_sza(n_days::I,
                            n_lats::I,
                            sza_arr::Array{FT},
                            scmap,
                            stitle,
                            file_name) where {FT<:AbstractFloat,I<:Int}
"""
function plot_day_lat_sza(d_arr::Array{I},
                                 l_arr::Array{FT},
                                 sza_arr::Array{FT},
                                 scmap,
                                 stitle,
                                 file_name) where {FT<:AbstractFloat,I<:Int}
  if scmap == "YlOrRd"
    cmap = :YlOrRd
    vmin, vmax = floor(min(sza_arr...)/10)*10, ceil(max(sza_arr...)/10)*10
  elseif scmap == "PRGn"
    cmap = :PRGn
    vmin, vmax = ceil(max(abs.(sza_arr)...)/10)*-10, ceil(max(abs.(sza_arr)...)/10)*10
  end
  
  p = contourf(d_arr, l_arr, sza_arr', c=cmap, clims=(vmin,vmax), title=stitle, 
    xlabel="Days since Jan 1", ylabel="Latitude", colorbar_title="Solar Zenith Angle [deg]")

  # subplot(121)
  # Fbar = sum(F_arr, dims=1)./size(F_arr,1)
  # plot(Fbar,l_arr,"k-", xlabel="Average ToA Insolation [W/m^2]")
  savefig(file_name)
end

# function main()
#   γ0 = 23.44
#   ϖ0 = 282.95
#   e0 = 0.017

#   days, lats, sza = calc_day_lat_sza(365, 180, γ0, ϖ0, e0)
#   title = "obliq=" * "$(γ0)" * ", perihelion=" * "$(ϖ0)" * ", ecc=" * "$(e0)"
#   plot_day_lat_sza(days, lats, sza, "YlOrRd", title, "example1.png")
# end

# main()
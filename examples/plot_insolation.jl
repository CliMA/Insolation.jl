using Insolation
using Insolation.SolarZenithAngle
using Plots

"""
    calc_day_lat_insolation(n_days::I,
                            n_lats::I,
                            γ::FT,
                            ϖ::FT,
                            e::FT) where {FT<:AbstractFloat,I<:Int}
"""
function calc_day_lat_insolation(n_days::I,
                                n_lats::I,
                                γ::FT,
                                ϖ::FT,
                                e::FT) where {FT<:AbstractFloat,I<:Int}
  d_arr = Array{I}(round.(collect(range(0, stop = 365, length = n_days))))
  l_arr = Array{I}(round.(collect(range(-90, stop = 90, length = n_lats))))
  F_arr = zeros(FT, n_days, n_lats)
  # loop over days
  for (i, day) in enumerate(d_arr)
    for (j, lat) in enumerate(l_arr)
      F_arr[i, j] = calc_point_insolation(day, lat, γ, ϖ, e)
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
                                 l_arr::Array{I},
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
    xlabel="Days since Jan 1", ylabel="Latitude", colorbar_title="ToA Insolation [W/m2]")

  # subplot(121)
  # Fbar = sum(F_arr, dims=1)./size(F_arr,1)
  # plot(Fbar,l_arr,"k-", xlabel="Average ToA Insolation [W/m^2]")
  savefig(file_name)
end
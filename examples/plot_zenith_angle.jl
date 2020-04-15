using Insolation
using Insolation.SolarZenithAngle

using Plots
# @static if Insolation.haspkg("Plots")
#   using Plots
#   const export_plots = true
# else
#   const export_plots = false
# end

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
  d_arr = Array{I}(round.(collect(range(0, stop = 365, length = n_day))))
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
                                 file_name) where {FT, I}
  if scmap == "jet" || scmap == "YlOrRd"
    vmin, vmax = 0, ceil(max(F_arr...)/100)*100
  else
    mm = ceil(max(abs.(F_arr)...)/10)*10
    vmin, vmax = -mm, mm
  end

  p1 = contourf(d_arr,l_arr,F_arr', title=stitle)
  plot!(p1, size=(400,400))
  p2 = contourf(d_arr,l_arr,F_arr', title=stitle)
  plot(p1,p2, layout=(1,2), xlabel="Days since Jan 1", ylabel="Latitude")
  plot!(p1, size=(400,400))
  savefig(file_name)

  # cbar = plt.colorbar()
  # cbar.set_label('ToA Insolation [W/m$^2$]', fontsize=f)

  # plt.subplot(gs[1])
  # Fbar = sum(F_arr, dims=1)./size(F_arr,1)
  # plt.plot(Fbar,l_arr,'k-')
  # plt.xlabel('Average ToA Insolation [W/m$^2$]', fontsize=f)
end
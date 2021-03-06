using Statistics

rtol = 1e-2

## Test globally averaged insolation ≈ TSI
ndays, nlats = [365, 181]
d_arr = Array{Int}(round.(collect(range(0, stop = 365, length = ndays))))
l_arr = collect(range(-90, stop = 90, length = nlats))
F_arr = zeros(ndays, nlats)

for (i, d) in enumerate(d_arr)
    for (j, lat) in enumerate(l_arr)
        datei = Dates.DateTime(2020,1,1) + Dates.Day(d)
        θ, dist = daily_zenith_angle(datei, lat, param_set)
        F_arr[i, j] = insolation(θ, dist, param_set)
    end
end

zonal_mean_insol = mean(F_arr, dims=1)
area_fac = abs.(cosd.(l_arr))
global_mean_insol = sum(zonal_mean_insol * area_fac) / sum(area_fac)
@test global_mean_insol ≈ tot_solar_irrad(param_set) / 4 rtol=rtol

## Test symmetry of insolation at equinox
date = Dates.DateTime(2021,3,20,9,34)
F_arr = zeros(nlats)
for (i, lat) in enumerate(l_arr)
    θ, dist = daily_zenith_angle(date, lat, param_set)
    F_arr[i] = insolation(θ, dist, param_set)
end
F_NH = sort(F_arr[l_arr .>= 0])
F_SH = sort(F_arr[l_arr .<= 0])
@test F_NH ≈ F_SH rtol=rtol

## Test invariance of zonal-mean insolation under rotation of ϖ
CLIMAParameters.Planet.lon_perihelion(::EarthParameterSet) = lon_perihelion_epoch(param_set) + π

ndays, nlats = [365, 181]
d_arr = Array{Int}(round.(collect(range(0, stop = 365, length = ndays))))
l_arr = collect(range(-90, stop = 90, length = nlats))
F_arr = zeros(ndays, nlats)

for (i, d) in enumerate(d_arr)
    for (j, lat) in enumerate(l_arr)
        datei = Dates.DateTime(2020,1,1) + Dates.Day(d)
        θ, dist = daily_zenith_angle(datei, lat, param_set)
        F_arr[i, j] = insolation(θ, dist, param_set)
    end
end

zonal_mean_insol_rotate = mean(F_arr, dims=1)
@test zonal_mean_insol_rotate ≈ zonal_mean_insol rtol=rtol
rtol = 0.1
atol = 1e-4

climlab_Insol = transpose(readdlm("climlab_daily_insolation.txt", ' ', FT, '\n'));

ndays, nlats = size(climlab_Insol)
d_arr = Array{Int}(round.(collect(range(0, stop = 365, length = ndays))))
l_arr = FT.(collect(range(-90, stop = 90, length = nlats)))
F_arr = zeros(ndays, nlats)

for (i, d) in enumerate(d_arr)
    for (j, lat) in enumerate(l_arr)
        datei = Dates.DateTime(2000,1,1) + Dates.Day(d)
        θ, dist = daily_zenith_angle(datei, lat, param_set)
        F_arr[i, j] = insolation(θ, dist, param_set)
    end
end

@test climlab_Insol ≈ F_arr rtol=rtol atol=atol

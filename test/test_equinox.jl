@testset "Equinox Date Tests" begin
    # Difference in NH and SH zenith angles at time x in given year
    function zdiff(x, year, od)
        date = xtomarchdate(x, year)
        Δt_years = Insolation.years_since_epoch(param_set, date)
        ϖ, γ, e = orbital_params(od, Δt_years)
        orb_params = (FT(ϖ), FT(γ), FT(e))
        result_s = daily_distance_zenith_angle(date, FT(-45), orb_params, param_set)
        result_n = daily_distance_zenith_angle(date, FT(45), orb_params, param_set)
        return result_n.daily_θ - result_s.daily_θ
    end

    # x is date relative to March 1, with 1.00 representing March 1 00:00
    function xtomarchdate(x, year)
        basedate = Dates.DateTime(year, 3, 1)
        deltat = Dates.Second(round((x - 1) * IP.day(param_set)))
        return basedate + deltat
    end

    od = Insolation.OrbitalDataSplines()
    days = zeros(length(1900:2100))
    for (i, year) in enumerate(1900:2100)
        f = (x -> zdiff(x, year, od))
        days[i] = find_zeros(f, 1.0, 30)[1]
    end

    # test mean is about March 21
    @test mean(days) ≈ 21 atol = 1

    # test decreasing
    @test mean(days[:100]) > mean(days[100:end])
end

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

    # The equinox date swings over a full day within this window (the leap cycle), so its
    # 1900-2100 mean is the meaningful quantity: March 20.70. The tolerance is set by the
    # Float32/Float64 spread of about 0.05 d, which dominates; plausible refinements of the
    # orbital elements move the mean by only a few thousandths of a day, so this stays a
    # real constraint on the calendar without tracking the parameter values.
    @test mean(days) ≈ 20.7 atol = 0.2

    # Equinox drifts to earlier dates within a century (a four-year leap cycle averages
    # 365.25 d, longer than the tropical year), reset by the century exceptions.
    @test mean(days[1:100]) > mean(days[100:end])
end

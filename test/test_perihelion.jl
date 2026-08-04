@testset "Perihelion Date Tests" begin
    # x is date relative to Jan 1, with 1.00 representing Jan 1 00:00
    function xtojandate(x, year)
        basedate = Dates.DateTime(year, 1, 1)
        deltat = Dates.Second(round((x - 1) * IP.day(param_set)))
        date = basedate + deltat
        return date
    end

    # Earth-Sun distance
    function edist(x, year, od)
        date = xtojandate(x, year)
        Δt_years = Insolation.years_since_epoch(param_set, date)
        ϖ, γ, e = orbital_params(od, Δt_years)
        orb_params = (FT(ϖ), FT(γ), FT(e))
        result = daily_distance_zenith_angle(date, FT(0), orb_params, param_set)
        # result.d is already in units of the semi-major axis
        return result.d
    end

    years = 1900:2100
    days = zeros(length(years))
    od = Insolation.OrbitalDataSplines()
    for (i, year) in enumerate(years)
        f = (x -> edist(x, year, od))
        res = optimize(f, 1.0, 30)
        days[i] = Optim.minimizer(res)[1]
    end

    # Mean date of perihelion over 1900-2100 is January 3.63. As for the equinox, the
    # tolerance is set by the Float32/Float64 spread (about 0.07 d) rather than by the
    # orbital elements, whose plausible refinements move the mean by ~0.003 d.
    @test mean(days) ≈ 3.6 atol = 0.2

    # Perihelion drifts to later dates: the Gregorian year is ~0.017 d shorter than the
    # anomalistic year, so perihelion arrives about 1.7 d later per century.
    @test mean(days[1:100]) < mean(days[100:end])
end

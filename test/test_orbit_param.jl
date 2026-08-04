od = Insolation.OrbitalDataSplines()

# The epoch parameters (supplied by ClimaParams) and the Laskar et al. (2004) tables are
# independent determinations of Earth's orbital elements at J2000, so evaluating the
# tables at t = 0 must reproduce the epoch parameters.
#
# The tolerances below are absolute and are derived from how fast each element actually
# changes, not from the current difference between the two sources. That keeps them stable
# when either source is refined, while still catching an element taken from the wrong
# epoch, given in the wrong unit, or belonging to a different body. Each tolerance is
# stated as the drift it corresponds to; the two sources currently agree far inside all
# three.
atol_ϖ = deg2rad(0.1)   # ϖ precesses ~1.7°/century, so this pins J2000 to within ~6 yr
atol_γ = deg2rad(0.02)  # γ drifts ~0.013°/century, so this pins J2000 to within ~1.5 cyr
atol_e = 1e-4           # e spans [0, 0.06] over the record; 0.6% of that range

@test mod(od.ϖ_spline(0.0), 2π) ≈ mod(IP.lon_perihelion_epoch(param_set), 2π) atol = atol_ϖ
@test od.γ_spline(0.0) ≈ IP.obliq_epoch(param_set) atol = atol_γ
@test od.e_spline(0.0) ≈ IP.eccentricity_epoch(param_set) atol = atol_e

ϖ0, γ0, e0 = orbital_params(od, 0.0)
@test mod(ϖ0, 2π) ≈ mod(IP.lon_perihelion_epoch(param_set), 2π) atol = atol_ϖ
@test γ0 ≈ IP.obliq_epoch(param_set) atol = atol_γ
@test e0 ≈ IP.eccentricity_epoch(param_set) atol = atol_e

# `year_anom` must be the anomalistic year, 365.2596 d. The competing definitions are the
# sidereal (365.2564 d), Julian (365.25 d), and tropical (365.2422 d) year; the nearest of
# them is 0.0032 d away, so the tolerance below separates them while remaining insensitive
# to refinements of the anomalistic year itself. Using the wrong one drifts a calendar date
# relative to perihelion by a full day every few thousand years of paleoclimate
# integration.
@test IP.year_anom(param_set) / IP.day(param_set) ≈ 365.2596 atol = 0.002

# Test broadcasting behavior
# Base.broadcastable(x::OrbitalDataSplines) = tuple(x) should make OrbitalDataSplines
# behave as a scalar in broadcasting operations
@testset "OrbitalDataSplines Broadcasting" begin
    # Test that OrbitalDataSplines is treated as a scalar in broadcasting
    @test Base.broadcastable(od) === tuple(od)

    # Test broadcasting with multiple time points
    times = [0.0, 1000.0, 5000.0]

    # This should work because od is treated as a scalar (broadcasts over times)
    results = orbital_params.(Ref(od), times)
    @test length(results) == 3

    # Each result should be iterable with 3 elements
    @test all(r -> length(collect(r)) == 3, results)

    # Verify the values are different for different times
    ϖ_values = [collect(r)[1] for r in results]
    @test !all(ϖ_values .≈ ϖ_values[1])  # Values should change with time

    # Test that we can collect and stack the results (common pattern in docs)
    y = hcat(collect.(orbital_params.(Ref(od), times))...)
    @test size(y) == (3, length(times))

    # Verify individual parameter extraction
    ϖ_arr, γ_arr, e_arr = y[1, :], y[2, :], y[3, :]
    @test length(ϖ_arr) == length(times)
    @test length(γ_arr) == length(times)
    @test length(e_arr) == length(times)

    # Check that all values are numeric (Real numbers)
    @test all(x -> x isa Real, ϖ_arr)
    @test all(x -> x isa Real, γ_arr)
    @test all(x -> x isa Real, e_arr)
end

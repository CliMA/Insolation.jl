rtol = 1e-2
od = Insolation.OrbitalDataSplines()
@test mod(od.ϖ_spline(0.0), 2π) ≈ mod(IP.lon_perihelion_epoch(param_set), 2π) rtol = rtol
@test od.γ_spline(0.0) ≈ IP.obliq_epoch(param_set) rtol = rtol
@test od.e_spline(0.0) ≈ IP.eccentricity_epoch(param_set) rtol = rtol

ϖ0, γ0, e0 = orbital_params(od, 0.0)
@test mod(ϖ0, 2π) ≈ mod(IP.lon_perihelion_epoch(param_set), 2π) rtol = rtol
@test γ0 ≈ IP.obliq_epoch(param_set) rtol = rtol
@test e0 ≈ IP.eccentricity_epoch(param_set) rtol = rtol

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

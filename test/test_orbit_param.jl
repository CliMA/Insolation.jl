rtol = 1e-2

@test mod(Insolation.ϖ_spline(0.0),2π) ≈ mod(IP.lon_perihelion_epoch(param_set),2π) rtol=rtol
@test Insolation.γ_spline(0.0) ≈ IP.obliq_epoch(param_set) rtol=rtol
@test Insolation.e_spline(0.0) ≈ IP.eccentricity_epoch(param_set) rtol=rtol

ϖ0, γ0, e0 = orbital_params(0.0)
@test mod(ϖ0,2π) ≈ mod(IP.lon_perihelion_epoch(param_set),2π) rtol=rtol
@test γ0 ≈ IP.obliq_epoch(param_set) rtol=rtol
@test e0 ≈ IP.eccentricity_epoch(param_set) rtol=rtol
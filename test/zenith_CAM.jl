# ϕ=latitude, θ=longitude, day=day in year since Jan 1
function CAM(ϕ, θ, day, year_since_1950)
    ϵ1 = deg2rad(23.320556)
    Aj = 0
    fj = 0
    δj = 0
    ϵ = ϵ1 + sum(Aj * cos(fj*year_since_1950 + δj))
    
    Mj = 0.01
    gj = 0
    βj = 0
    RHScos = sum(Mj*cos(gj*year_since_1950 + βj))
    RHSsin = sum(Mj*sin(gj*year_since_1950 + βj))
    e = sqrt(RHScos^2 + RHSsin^2)
    Π = atan(RHScos/RHSsin)
    
    ψ1 = deg2rad(50.439273 / 360 / 360)
    ξ = deg2rad(3.392506)
    Fj = 0
    fj_ = 0
    δj_ = 0
    ψ = ψ1*year_since_1950 + ξ + sum(Fj*sin(fj_*year_since_1950 + δj_))
    ϖ = Π + ψ
    
    β = sqrt(1-e^2)
    λm0 = (e + e^3/4)*(1+β)*sin(ϖ) - (e^2/2)*(1/2 + β)*sin(2ϖ) + (e^3/4)*(1/3 + β)*sin(3ϖ)
    dve = 80.5
    λm = λm0 + 2π*(day-dve) / 365
    λ = λm + (2e - e^3/4)*sin(λm - ϖ) + (5e^2/4)*sin(2*(λm - ϖ)) + (13e^3/12)*sin(3*(λm - ϖ))

    δ = asin(sin(ϵ)*sin(λ))
    H = 2π * (day + θ/2π)
    cosμ = sin(ϕ)*sin(δ) - cos(ϕ)*cos(δ)*cos(H)
    
    ρ = (1 - e^2) / (1 + e*cos(λ-ϖ))
    μ = rad2deg(acos(cosμ))

    # println(ϵ)
    # println(e)
    # println(Π)
    # println(ϖ)
    # println(λ)
    println(μ, "\t", ρ)

    return μ, ρ
end

X = 10
CAM(0, 0, X + 0.0, 75);
CAM(0, 0, X + 0.3, 75);
CAM(0, 0, X + 0.5, 75);
CAM(0, 0, X + 0.7, 75);
CAM(0, 0, X + 1.0, 75);

X = 100
CAM(0, 0, X + 0.0, 75);
CAM(0, 0, X + 0.3, 75);
CAM(0, 0, X + 0.5, 75);
CAM(0, 0, X + 0.7, 75);
CAM(0, 0, X + 1.0, 75);

X = 200
CAM(0, 0, X + 0.0, 75);
CAM(0, 0, X + 0.3, 75);
CAM(0, 0, X + 0.5, 75);
CAM(0, 0, X + 0.7, 75);
CAM(0, 0, X + 1.0, 75);
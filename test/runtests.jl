push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Test

using Insolation 
using Dates
using Statistics
using Roots
using Optim

using CLIMAParameters
using CLIMAParameters.Planet: orbit_semimaj, tot_solar_irrad
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

ϖ_spline, γ_spline, e_spline = orbital_params_spline();
CLIMAParameters.Planet.lon_perihelion_spline(::EarthParameterSet) = ϖ_spline;
CLIMAParameters.Planet.obliq_spline(::EarthParameterSet) = γ_spline;
CLIMAParameters.Planet.eccentricity_spline(::EarthParameterSet) = e_spline;

@testset "Types" begin
    include("test_types.jl")
end

FT = Float64
@testset "Zenith Angle" begin
    include("test_zenith_angle.jl")
end
@testset "Insolation" begin
    include("test_insolation.jl")
end
@testset "Orbits" begin
    include("test_perihelion.jl")
    include("test_equinox.jl")
end

module Insolation

using Dates, DelimitedFiles, Interpolations

include("Parameters.jl")
import .Parameters as IP
const AIP = IP.AbstractInsolationParams

export orbital_params

mock_t_range = 1.0:1.0:10.0;
mock_x_data = rand(10);
const e_spline_etp = Ref(CubicSplineInterpolation(mock_t_range, mock_x_data, extrapolation_bc = NaN))
const γ_spline_etp = Ref(CubicSplineInterpolation(mock_t_range, mock_x_data, extrapolation_bc = NaN))
const ϖ_spline_etp = Ref(CubicSplineInterpolation(mock_t_range, mock_x_data, extrapolation_bc = NaN))
function __init__()
    datapath = joinpath(@__DIR__, "../src/data/", "INSOL.LA2004.BTL.csv");
    x, _ = readdlm(datapath, ',', Float64, header=true);
    t_range = x[1,1]*1e3 : 1e3 : x[end,1]*1e3; # array of every 1 kyr to range of years
    e_spline_etp[] = CubicSplineInterpolation(t_range, x[:,2], extrapolation_bc = NaN);
    γ_spline_etp[] = CubicSplineInterpolation(t_range, x[:,3], extrapolation_bc = NaN);
    ϖ_spline_etp[] = CubicSplineInterpolation(t_range, x[:,4], extrapolation_bc = NaN);
end
function e_spline(args...)
    return e_spline_etp[](args...)
end
function γ_spline(args...)
    return γ_spline_etp[](args...)
end
function ϖ_spline(args...)
    return ϖ_spline_etp[](args...)
end

"""
    orbital_params(dt::FT) where {FT <: Real}

This function returns the orbital parameters (ϖ, γ, e) at dt (years) since J2000 epoch.
Data are read from file and interpolation functions are created in __init__() method.
The functions are stored as global variables that are used inside Insolation.jl.
The parameters vary due to Milankovitch cycles. 
Orbital parameters from the Laskar 2004 paper are in the "src/data/INSOL.LA2004.BTL.csv" file.
"""
function orbital_params(dt::FT) where {FT <: Real}
    return ϖ_spline(dt), γ_spline(dt), e_spline(dt)
end

include("ZenithAngleCalc.jl")
include("InsolationCalc.jl")

end # module
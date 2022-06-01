module Insolation

using Dates, DelimitedFiles, Interpolations
using CLIMAParameters
using CLIMAParameters.Planet
using CLIMAParameters: AbstractParameterSet
const APS = AbstractParameterSet

export orbital_params

"""
This function reads orbital parameters from the Laskar 2004 data set and creates 
interpolation functions for the parameters as a function of the years since the J2000 epoch.
The functions are stored as global variables that are used inside of Insolation.jl.
The parameters vary due to Milankovitch cycles. 
Data from this paper are in the "src/data/INSOL.LA2004.BTL.csv" file.
"""
function __init__()
    datapath = joinpath(@__DIR__, "../src/data/", "INSOL.LA2004.BTL.csv");
    x, _ = readdlm(datapath,',',header=true);
    t_step, t_units = 1, 1e3;
    t_range = x[1,1]*t_units : t_step*t_units : x[end,1]*t_units;
    
    global e_spline = CubicSplineInterpolation(t_range, x[:,2], extrapolation_bc = NaN);
    global γ_spline = CubicSplineInterpolation(t_range, x[:,3], extrapolation_bc = NaN);
    global ϖ_spline = CubicSplineInterpolation(t_range, x[:,4], extrapolation_bc = NaN);
end

"""
    orbital_params(dt::FT) where {FT <: Real}

This function returns the orbital parameters (ϖ, γ, e) at dt (years) since J2000 epoch.
Data are read from file and interpolation function are created in __init__() of Insolation.jl
"""
function orbital_params(dt::FT) where {FT <: Real}
    return ϖ_spline(dt), γ_spline(dt), e_spline(dt)
end

include("ZenithAngleCalc.jl")
include("InsolationCalc.jl")

end # module
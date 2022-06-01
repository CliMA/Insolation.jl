using DelimitedFiles, Interpolations

export orbital_params_spline

"""
    orbital_params_spline()

This function returns functions for the orbital parameters (ϖ, γ, e) 
which are functions of the years post J2000 epoch.
The parameters vary due to Milankovitch cycles. 
The dominant 10 frequencies of these cycles are used based on a
Fourier analysis of the parameters as calculated in the
Lasker et al. (2004) paper.
Data from this paper are in the "src/data/INSOL.LA2004.BTL.csv" file.
"""
function orbital_params_spline()
    datapath = joinpath(@__DIR__, "../src/data/", "INSOL.LA2004.BTL.csv");
    x, _ = readdlm(datapath,',',header=true);
    t_step, t_units = 1, 1e3;
    t_range = x[1,1]*t_units : t_step*t_units : x[end,1]*t_units;
    
    e_spline = CubicSplineInterpolation(t_range, x[:,2], extrapolation_bc = NaN);
    γ_spline = CubicSplineInterpolation(t_range, x[:,3], extrapolation_bc = NaN);
    ϖ_spline = CubicSplineInterpolation(t_range, x[:,4], extrapolation_bc = NaN);
    return [ϖ_spline, γ_spline, e_spline]
end

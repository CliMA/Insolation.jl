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
    x,header = readdlm(datapath,',',header=true);
    t_step, t_units = 1, 1e3;
    t_range = x[1,1]*t_units : t_step*t_units : x[end,1]*t_units;
    
    e_spline = CubicSplineInterpolation(t_range, x[:,2]);
    γ_spline = CubicSplineInterpolation(t_range, x[:,3]);
    ϖ_spline = ConstantInterpolation(t_range, x[:,4]);
    ## this needs to be improved but I'm having some trouble getting a higher order 
    ## interpolation to not return something outside the valid range
    ## tried doing the interpolation on sinϖ, that might be a way forward?
    return [ϖ_spline, γ_spline, e_spline]
end

module Insolation

using Dates, DelimitedFiles, Interpolations
import ArtifactWrappers as AW

include("Parameters.jl")
import .Parameters as IP
const AIP = IP.AbstractInsolationParams

export orbital_params

function orbital_parameters_dataset_path()
    era_dataset = AW.ArtifactWrapper(
        @__DIR__,
        "era-global",
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/3y02smnlxhgwednm3eho7lve2xhq1n7r.csv",
            filename = "INSOL.LA2004.BTL.csv",
        ),],
    )
    return AW.get_data_folder(era_dataset)
end

"""
    OrbitalData

The parameters vary due to Milankovitch cycles. 

Orbital parameters from the Laskar 2004 paper are
lazily downloaded from Caltech Box to the
`orbital_parameters_dataset_path()` path.
"""
struct OrbitalData{E, G, O}
    e_spline_etp::E
    γ_spline_etp::G
    ϖ_spline_etp::O
    function OrbitalData()
        datapath = joinpath(orbital_parameters_dataset_path(), "INSOL.LA2004.BTL.csv");
        x, _ = readdlm(datapath, ',', Float64, header=true);
        t_range = x[1,1]*1e3 : 1e3 : x[end,1]*1e3; # array of every 1 kyr to range of years
        e_spline_etp = CubicSplineInterpolation(t_range, x[:,2], extrapolation_bc = NaN);
        γ_spline_etp = CubicSplineInterpolation(t_range, x[:,3], extrapolation_bc = NaN);
        ϖ_spline_etp = CubicSplineInterpolation(t_range, x[:,4], extrapolation_bc = NaN);

        E = typeof(e_spline_etp)
        G = typeof(γ_spline_etp)
        O = typeof(ϖ_spline_etp)
        return new{E,G,O}(e_spline_etp,γ_spline_etp,ϖ_spline_etp)
    end
end

Base.broadcastable(x::OrbitalData) = tuple(x)

e_spline(od, args...) = od.e_spline_etp(args...)
γ_spline(od, args...) = od.γ_spline_etp(args...)
ϖ_spline(od, args...) = od.ϖ_spline_etp(args...)

"""
    orbital_params(od::OrbitalData, dt::FT) where {FT <: Real}

Parameters are interpolated from the values given in the
Laskar 2004 dataset using a cubic spline interpolation.

See [`OrbitalData`](@ref).
"""
function orbital_params(od::OrbitalData, dt::FT) where {FT <: Real}
    return ϖ_spline(od, dt), γ_spline(od, dt), e_spline(od, dt)
end

include("ZenithAngleCalc.jl")
include("InsolationCalc.jl")

end # module
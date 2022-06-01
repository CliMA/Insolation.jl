module Insolation

using Dates
using CLIMAParameters: AbstractParameterSet
using CLIMAParameters
using CLIMAParameters.Planet
const APS = AbstractParameterSet

include("ZenithAngleCalc.jl")
include("InsolationCalc.jl")
include("OrbitalParams.jl")

end # module
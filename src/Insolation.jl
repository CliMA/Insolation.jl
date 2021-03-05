module Insolation

using Dates

using CLIMAParameters: AbstractParameterSet
using CLIMAParameters
using CLIMAParameters.Planet
const APS = AbstractParameterSet
Base.broadcastable(param_set::APS) = Ref(param_set)

include("ZenithAngleCalc.jl")
include("InsolationCalc.jl")

end # module
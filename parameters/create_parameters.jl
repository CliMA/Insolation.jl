import CLIMAParameters as CP
import Insolation.Parameters as IP

function create_insolation_parameters(FT, overrides::NamedTuple = NamedTuple())
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    aliases = string.(fieldnames(IP.InsolationParameters))
    pairs = CP.get_parameter_values!(toml_dict, aliases, "Insolation")
    params = merge((;pairs...), overrides) # overrides
    param_set = IP.InsolationParameters{FT}(; params...)
    return param_set
end

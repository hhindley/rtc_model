using DifferentialEquations, DataFrames, ModelingToolkit, OrderedCollections

PATH = "YOURPATH"

include("$PATH/funcs.jl");
include("$PATH/params.jl");
include("$PATH/models/rtc_model.jl");
include("$PATH/models/rtc_inhibition_model.jl");
include("$PATH/models/rtc_trna_model.jl");
include("$PATH/models/trna_inhib_models.jl");

# rtc model
solu_rtc = sol(rtc_model, init_rtc, tspan, params_rtc)

# rtc model with damage 
params_dam = deepcopy(params_rtc)
params_dam[kdam] = 0.1
solu_rtc_dam = sol(rtc_model, ssvals_rtc, tspan, params_dam)

# rtcb inhibition model 
solu_rtc_inhib = sol(rtcb_inhib_model, init_inhib_rtcb, tspan, params_inhib)

# trna model 
solu_trna = sol(rtc_trna_model, init_trna, tspan, params_trna)

# trna model with damage 
params_dam_trna = deepcopy(params_trna)
params_dam_trna[kdam] = 0.1
solu_trna_dam = sol(rtc_trna_model, ssvals_trna, tspan, params_dam_trna)

# rtcb trna inhibition model
solu_trna_inhib = sol(rtcb_trna_inhib_model, init_trna_inhib_rtcb, tspan, params_trna_inhib)


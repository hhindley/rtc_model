using DifferentialEquations, DataFrames, ModelingToolkit, OrderedCollections

include(joinpath(homedir(), "rtc_model/funcs.jl"));
include(joinpath(homedir(), "rtc_model/rtc_params.jl"));
include(joinpath(homedir(), "rtc_model/trna_params.jl"));

# rRNA models 
include(joinpath(homedir(), "rtc_model/models/rtc_orig.jl"));
include(joinpath(homedir(), "rtc_model/models/rtc_inhibition_model.jl"));
include(joinpath(homedir(), "rtc_model/models/rtc_orig_cooperativity.jl"));
include(joinpath(homedir(), "rtc_model/models/rtc_without_rtca.jl"));
include(joinpath(homedir(), "rtc_model/models/lam_kin.jl"));

# tRNA models
include(joinpath(homedir(), "rtc_model/models/rtc_trna_model.jl"));
include(joinpath(homedir(), "rtc_model/models/trna_inhib_models.jl"));
include(joinpath(homedir(), "rtc_model/models/rtc_trna_cooperativity.jl"));
include(joinpath(homedir(), "rtc_model/models/rtc_trna_model_without_rtca.jl"));

# below are some examples of using DifferentialEquations.jl to solve models, the same code structure can be 
# used for the other models loaded above.

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


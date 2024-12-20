using DifferentialEquations, DataFrames, ModelingToolkit, OrderedCollections, LinearAlgebra

include(joinpath(homedir(), "rtc_model/funcs.jl")) 
include(joinpath(homedir(), "rtc_model/rtc_params.jl"))
include(joinpath(homedir(), "rtc_model/trna_params.jl"))

# rRNA models 
include(joinpath(homedir(), "rtc_model/models/rtc_orig.jl"))
include(joinpath(homedir(), "rtc_model/models/rtc_inhibition_model.jl"))
include(joinpath(homedir(), "rtc_model/models/rtc_orig_cooperativity.jl"))
include(joinpath(homedir(), "rtc_model/models/rtc_without_rtca.jl"))
include(joinpath(homedir(), "rtc_model/models/lam_kin.jl"))

# tRNA models
include(joinpath(homedir(), "rtc_model/models/rtc_trna_model.jl"))
include(joinpath(homedir(), "rtc_model/models/trna_inhib_models.jl"))
include(joinpath(homedir(), "rtc_model/models/rtc_trna_cooperativity.jl"))
include(joinpath(homedir(), "rtc_model/models/rtc_trna_model_without_rtca.jl"))

# this must be loaded after the models set up
using BifurcationKit

# below are two examples of using BifurcationKit.jl to find bifurcation points, the same code structure can be 
# used for the other models loaded above.

# rtc model 
br = get_br(rtc_model , ssvals_rtc, params_rtc, kdam_max=0.8)

# trna model 
br_trna = get_br(rtc_trna_model, ssvals_trna, params_trna, kdam_max=0.8)


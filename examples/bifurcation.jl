using DifferentialEquations, DataFrames, ModelingToolkit, OrderedCollections, LinearAlgebra

include(joinpath(homedir(), "rtc_model/src/funcs.jl")) 
include(joinpath(homedir(), "rtc_model/src/parameters/rtc_params.jl"))
include(joinpath(homedir(), "rtc_model/src/parameters/trna_params.jl"))

# rRNA models 
include(joinpath(homedir(), "rtc_model/src/models/rtc_orig.jl"))
include(joinpath(homedir(), "rtc_model/src/models/rtc_inhibition_model.jl"))
include(joinpath(homedir(), "rtc_model/src/models/rtc_orig_cooperativity.jl"))
include(joinpath(homedir(), "rtc_model/src/models/rtc_without_rtca.jl"))
include(joinpath(homedir(), "rtc_model/src/models/lam_kin.jl"))

# tRNA models
include(joinpath(homedir(), "rtc_model/src/models/rtc_trna_model.jl"))
include(joinpath(homedir(), "rtc_model/src/models/trna_inhib_models.jl"))
include(joinpath(homedir(), "rtc_model/src/models/rtc_trna_cooperativity.jl"))
include(joinpath(homedir(), "rtc_model/src/models/rtc_trna_model_without_rtca.jl"))

# this must be loaded after the models set up
using BifurcationKit

# below are two examples of using BifurcationKit.jl to find bifurcation points, the same code structure can be 
# used for the other models loaded above.

# rtc model 
br = get_br(rtc_model , ssvals_rtc, params_rtc, kdam_max=0.8)

# trna model 
br_trna = get_br(rtc_trna_model, ssvals_trna, params_trna, kdam_max=0.8)

# example of plotting a bifurcation diagram
using Plots
df = create_br_df(br)
plot(df.kdam, df.rtca)
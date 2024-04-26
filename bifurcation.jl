using DifferentialEquations, DataFrames, ModelingToolkit, OrderedCollections, LinearAlgebra

PATH = "YOURPATH"

include("$PATH/funcs.jl");
include("$PATH/params.jl");
include("$PATH/models/rtc_model.jl");
include("$PATH/models/rtc_inhibition_model.jl");
include("$PATH/models/rtc_trna_model.jl");
include("$PATH/models/trna_inhib_models.jl");

# rtc model 
br = get_br(rtc_model, ssvals_rtc, params_rtc, 1.5)

# trna model 
br = get_br(rtc_trna_model, ssvals_trna, params_trna, 20.)


# A model of Rtc-regulated RNA repair reveals heterogeneous responses to ribosome-targeting antibiotics and suggests molecular targets to potentiate their effects

## Provided here is the code defining and solving the models used in the paper.

Provided you download these files into your machines home directory, the files should run with the below dependencies installed. run_model.jl and bifurcation.jl are demo files.

### Software dependencies
- Analyses performed on operating system Ubuntu 22.04.2 LTS or MacOS 14.7.2 
- Julia v1.11
- ModelingToolkit.jl (v9.54.0) 
- DifferentialEquations.jl (v7.15.0)
- BifurcationKit.jl (v0.4.4)

### Installation guide
Refer to https://julialang.org/downloads/ for instructions on how to download julia and packages.

### models folder
- includes all model versions.
  
### funcs.jl 
- includes all functions needed that are used within other files.
  
### params.jl 
- includes all parameter values as used to get results.

### run_model.jl 
- provides all information needed to solve each model with and without damage.
  
### bifurcation.jl 
- provides all information needed to obtain the bistability result in the Rtc and tRNA models.
- when plotting these results for each species the outputs should match what is shown in paper figure 3a.
- file should not take more than a few minutes to run on a standard machine.


# Heterogeneity in responses to ribosome-targeting antibiotics mediated by bacterial RNA repair

## Installation guide
Refer to https://julialang.org/downloads/ for instructions on how to download Julia.

To get the code, clone the repository:
```sh
git clone https://github.com/hhindley/rtc_model.git
cd rtc_model
```

Provided you download these files into your machine's home directory, the files should run with the below dependencies installed. 

## Julia Environment

This project includes a `Project.toml` file to specify required packages. The `Manifest.toml` is excluded from version control because it is machine-specific.  
To set up the environment and generate your own `Manifest.toml`, run the following in the Julia REPL from the project directory:

```julia
] activate .
] instantiate
```

This will install all required packages and create a local `Manifest.toml` for reproducibility.

### Software dependencies
- Analyses performed on operating system Ubuntu 22.04.2 LTS or MacOS 14.7.2 
- Julia v1.11
- ModelingToolkit.jl (v9.54.0) 
- DifferentialEquations.jl (v7.15.0)
- BifurcationKit.jl (v0.4.4)

These dependencies will be directly installed if you follow the above instructions to setup the julia environment. 

## Directory Structure

```
rtc_model/
├── Project.toml
├── README.md
├── examples/
│   ├── bifurcation.jl
│   └── run_model.jl
├── src/
│   ├── funcs.jl
│   ├── models/
│   │   ├── lam_kin.jl
│   │   ├── rtc_inhibition_model.jl
│   │   ├── rtc_orig.jl
│   │   ├── rtc_orig_cooperativity.jl
│   │   ├── rtc_trna_cooperativity.jl
│   │   ├── rtc_trna_model.jl
│   │   ├── rtc_trna_model_without_rtca.jl
│   │   ├── rtc_without_rtca.jl
│   │   └── trna_inhib_models.jl
│   ├── parameters/
│   │   ├── rtc_params.jl
│   │   └── trna_params.jl
```

### examples/
- `run_model.jl`: Provides all information needed to solve each model with and without damage.
- `bifurcation.jl`: Provides all information needed to obtain the bistability result in the Rtc and tRNA models. When plotting these results for each species, the outputs should match what is shown in paper figure 3a. Should not take more than a few minutes to run on a standard machine.

### src/
- `funcs.jl`: Includes all functions used within other files.
- `models/`: Contains all model versions.
- `parameters/`: Contains all parameter values as used to get results.

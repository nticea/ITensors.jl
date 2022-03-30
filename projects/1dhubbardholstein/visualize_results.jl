using Pkg
Pkg.activate(joinpath(@__DIR__,"../.."))
include(joinpath(@__DIR__,"model.jl"))
include(joinpath(@__DIR__,"plotting.jl"))
include(joinpath(@__DIR__,"utilities.jl"))
using ITensors
using ITensors.HDF5

# Load the results in 
## NOTE: NEED TO CHANGE UTILS BACK TO TAKE IN PHONON ARGUMENTS TO PARAMS ## 
loadpath = "/Users/nicole/sherlock/ITensors.jl/projects/1dhubbardholstein/outputs/2022-03-22_10:18:27.h5"
params, hubbholst, dmrg_results, equilibrium_corr, tebd_results = load_structs(loadpath)

# Visualize things we want
plot_equilibrium_correlations(dmrg_results,"spin",hubbholst)
plot_entropy(tebd_results)
plot_correlation_function(tebd_results)
plot_spectral_function(tebd_results,params)




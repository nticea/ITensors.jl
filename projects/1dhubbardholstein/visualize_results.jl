using Pkg
Pkg.activate(joinpath(@__DIR__,"../.."))
include(joinpath(@__DIR__,"model.jl"))
include(joinpath(@__DIR__,"plotting.jl"))
include(joinpath(@__DIR__,"utilities.jl"))
using ITensors
using ITensors.HDF5

# Load the results in 
fname = "2022-03-20_10:25:17" # format yyyy-mm-dd_HH:MM:SS
loadpath = joinpath(@__DIR__,"outputs",fname*".h5")
params, hubbholst, dmrg_results, equilibrium_corr, tebd_results = load_structs(loadpath)

# Visualize things we want
plot_equilibrium_correlations(dmrg_results,"spin",hubbholst)
plot_entropy(tebd_results)
plot_correlation_function(tebd_results)
plot_spectral_function(tebd_results,params)




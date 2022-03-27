## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__,"../.."))
using Dates
include(joinpath(@__DIR__,"model.jl"))
#include(joinpath(@__DIR__,"plotting.jl"))
include(joinpath(@__DIR__,"utilities.jl"))

## SAVING INFO ##
DO_SAVE = true
INTERIM_SAVE = true

## PARAMETERS ## 

# Model 
N = 80
t = 1 
U = 8
ω = 0*t 
g0 = 0*t 
g1 = 0*g0 
doping = 0
max_num_phonons = 1 ## TODO: incorporate this! ##

# Simulation 
T = 30
τ = 0.05
DMRG_numsweeps = 80
DMRG_maxdim = 800
TEBD_maxdim = 800
TEBD_cutoff = 1E-10
DMRG_cutoff = 1E-10

# Saveout info 
fname_out = Dates.format(now(), "no_doping_1_phonon_80")
save_path = joinpath(@__DIR__,"outputs",fname_out*".h5")

# Specify operators of interest
A_t0 = "Cup"
A_t = "Cdagup"

## CODE ## 

# Initialize 
println("Initializing...")
params = parameters(N=N, t=t, U=U, ω=ω, g0=g0, g1=g1, doping=doping, 
                    DMRG_numsweeps=DMRG_numsweeps,
                    DMRG_maxdim=DMRG_maxdim, DMRG_cutoff=DMRG_cutoff,
                    T=T, τ=τ, TEBD_cutoff=TEBD_cutoff)
hubbholst = HubbardHolsteinModel(params)
if DO_SAVE
    save_structs(params, save_path)
end

# Run DMRG
println("Finding ground state...")
dmrg_results = run_DMRG(hubbholst, params, alg="divide_and_conquer")
if DO_SAVE
    save_structs(dmrg_results, save_path)
end

# Equilibrium correlations
println("Computing equilibrium correlations...")
eq_corr = compute_all_equilibrium_correlations(dmrg_results, hubbholst)
if DO_SAVE
    save_structs(eq_corr, save_path)
end

# Compute correlation functions 
println("Computing correlation functions...")
tebd_results = compute_correlations(dmrg_results, A_t0, A_t, hubbholst, params,
                                    interim_save=INTERIM_SAVE, savepath=save_path)
if DO_SAVE
    save_structs(tebd_results, save_path)
end
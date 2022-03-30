## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__,"../.."))
using Dates
include(joinpath(@__DIR__,"model.jl"))
include(joinpath(@__DIR__,"utilities.jl"))

## SAVING INFO ##

DO_SAVE = false
INTERIM_SAVE = false

## PARAMETERS ## 

# Model 
N = 8
t = 1 
U = 8
ω = 0.1*t 
g0 = 0.5*t 
g1 = 0.1*t 
doping = 0
max_phonons = 2 

# save path
date_stamp = Dates.format(now(), "HH:MM:SS") 
param_stamp = "_$(N)N_$(t)t_$(U)U_$(ω)ω_$(g0)g0_$(g1)g1_$(doping)doping_$(max_phonons)phonons"
save_path = joinpath(@__DIR__,"outputs",date_stamp*param_stamp*".h5")

# Simulation 
T = 10
τ = 0.01
DMRG_numsweeps = 80
DMRG_maxdim = 800
TEBD_maxdim = 800
TEBD_cutoff = 1E-10
DMRG_cutoff = 1E-10

# Specify spectral function operators 
A_t0 = "Cup"
A_t = "Cdagup"

## CODE ## 

# Initialize 
println("Initializing...")
params = parameters(N=N, t=t, U=U, ω=ω, g0=g0, g1=g1, doping=doping, 
                    max_phonons=max_phonons,
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
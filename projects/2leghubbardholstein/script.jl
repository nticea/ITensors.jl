## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__,"../.."))
using Dates
include(joinpath(@__DIR__,"model.jl"))
include(joinpath(@__DIR__,"utilities.jl"))

## SAVING INFO ##
DO_SAVE = true
INTERIM_SAVE = true

## PARAMETERS ## 

# Model 
Nx = 20
Ny = 2
t = 1 
U = 8
ω = 0.5*t 
g0 = 0.1*t 
g1 = 0.05*t 
λ = 0*t
doping = 0
max_phonons = 5

# Save 
date_stamp = Dates.format(now(), "y-m-d_HH:MM:SS") 
param_stamp = "_$(Nx)Nx_$(Ny)Ny_$(t)t_$(U)U_$(ω)ω_$(g0)g0_$(g1)g1_$(λ)λ_$(doping)doping_$(max_phonons)phonons"
save_path = joinpath(@__DIR__,"outputs",date_stamp*param_stamp*".h5")

# Simulation 
DMRG_numsweeps = 80
DMRG_maxdim = 2000
DMRG_cutoff = 1E-10
DMRG_LBO = true
T = 1#80
τ = 0.05
TEBD_maxdim = 1000
TEBD_cutoff = 1E-10
TEBD_LBO = true

# Specify spectral function operators 
A_t0 = "Cup"
A_t = "Cdagup"

## CODE ## 

# Initialize 
println("Initializing...")
params = parameters(Nx=Nx, Ny=Ny, t=t, U=U, ω=ω, g0=g0, g1=g1, λ=λ, doping=doping, 
                    max_phonons=max_phonons,
                    DMRG_numsweeps=DMRG_numsweeps,
                    DMRG_maxdim=DMRG_maxdim, DMRG_cutoff=DMRG_cutoff,DMRG_LBO=DMRG_LBO,
                    T=T, τ=τ, TEBD_cutoff=TEBD_cutoff,TEBD_LBO=TEBD_LBO)
hubbholst = HubbardHolsteinModel(params)
if DO_SAVE
    save_structs(params, save_path)
    save_structs(hubbholst, save_path)
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
# println("Computing correlation functions...")
# tebd_results = compute_correlations(dmrg_results, A_t0, A_t, hubbholst, params,
#                                     interim_save=INTERIM_SAVE, savepath=save_path)
# if DO_SAVE
#     save_structs(tebd_results, save_path)
# end
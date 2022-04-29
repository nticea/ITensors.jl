## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__,"../.."))
include(joinpath(@__DIR__,"model.jl"))
include(joinpath(@__DIR__,"utilities.jl"))

## PARAMETERS ## 

# Model 
N = 8
t = 1 
U = 8
ω = 0.2*t 
g0 = 0.2*t 
g1 = 0.15*t 
λ = 0.05*t
doping = 0

# Simulation 
T = 40
τ = 0.05
DMRG_numsweeps = 20
DMRG_maxdim = 800
TEBD_maxdim = 800
TEBD_cutoff = 1E-10
DMRG_cutoff = 1E-10

# Number of phonons to try 
max_phonons_list = [1, 2, 3, 4, 5, 7]
gs_energy_λ = []
gs_energy_no_λ = []
phonon_density_λ = []
phonon_density_no_λ = []
ΔEs = []

for max_phonons in max_phonons_list
    
    # Initialize 
    println("Initializing with ", max_phonons, " phonons...")
    params = parameters(N=N, t=t, U=U, ω=ω, g0=g0, g1=g1, λ=λ, doping=doping, 
                        max_phonons=max_phonons,
                        DMRG_numsweeps=DMRG_numsweeps,
                        DMRG_maxdim=DMRG_maxdim, DMRG_cutoff=DMRG_cutoff,
                        T=T, τ=τ, TEBD_cutoff=TEBD_cutoff)
    hubbholst = HubbardHolsteinModel(params)

    # Run DMRG
    println("Finding ground state for λ=",λ)
    dmrg_results = run_DMRG(hubbholst, params, alg="divide_and_conquer")
    
    # compensate for offset that comes from having a nonzero λ
    E = dmrg_results.ground_state_energy
    nb = compute_phonon_number(dmrg_results.ground_state)
    ΔE = λ*sum(nb.^2)

    push!(ΔEs, ΔE)
    push!(gs_energy_λ, E)
    push!(phonon_density_λ, nb)

    # Initialize 
    println("Initializing with 0 λ")
    params = parameters(N=N, t=t, U=U, ω=ω, g0=g0, g1=g1, λ=0, doping=doping, 
                        max_phonons=max_phonons,
                        DMRG_numsweeps=DMRG_numsweeps,
                        DMRG_maxdim=DMRG_maxdim, DMRG_cutoff=DMRG_cutoff,
                        T=T, τ=τ, TEBD_cutoff=TEBD_cutoff)
    hubbholst = HubbardHolsteinModel(params)

    # Run DMRG
    println("Finding ground state for λ=0")
    dmrg_results = run_DMRG(hubbholst, params, alg="divide_and_conquer")
    
    # compensate for offset that comes from having a nonzero λ
    E = dmrg_results.ground_state_energy
    nb = compute_phonon_number(dmrg_results.ground_state)
    ΔE = λ*sum(nb.^2)

    push!(ΔEs, ΔE)
    push!(gs_energy_no_λ, E)
    push!(phonon_density_no_λ, nb)
end

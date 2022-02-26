using Pkg; 
Pkg.activate("././.")
using Revise
using ITensors

## FUNCTIONS ##

struct Parameters
    # Model parameters
    N::Int
    t::Real
    U::Real
    ω::Real
    g0::Real
    g1::Real
    doping::Real

    # DMRG parameters
    sweeps::Sweeps

    # TEBD parameters
    mid::Int
    T::Int
    τ::Real
    TEBD_cutoff
    TEBD_maxdim
end

struct HubbardHolsteinModel
    # Lattice
    sites 
    # Hamiltonian (in two forms)
    mpo::MPO
    gates 
    ampo
end

struct DMRGResults
    ground_state
    ground_state_energy
    entropy 
end

struct TEBDResults
    entropy
    self_overlap
    phonon_flux
    corrs 
end

function parameters(;N::Int, t::Real, U::Real=nothing, ω::Real=nothing, 
                    g0::Real=nothing, g1::Real=nothing, doping::Real=0, 
                    DMRG_numsweeps::Int=20, DMRG_noise=nothing, DMRG_maxdim=nothing, DMRG_cutoff=nothing,
                    T::Int=25, τ::Real=0.1, TEBD_cutoff=1E-14, TEBD_maxdim=400)
    
    if isnothing(U)
        U = 8*t
    end
    if isnothing(ω)
        ω = 0.5*t
    end
    if isnothing(g0)
        g0 = 0.05*t
    end
    if isnothing(g1)
        g1 = 0*g0
    end

    sweeps = Sweeps(DMRG_numsweeps)
    if isnothing(DMRG_noise)
        DMRG_noise = [1E-6,1E-6,1E-8,0]
    end
    if isnothing(DMRG_maxdim)
        DMRG_maxdim = [20,40,100,200,400]
    end

    if isnothing(DMRG_cutoff)
        DMRG_cutoff = 1E-14
    end
    setnoise!(sweeps, DMRG_noise...) # Very important to use noise for this model
    setmaxdim!(sweeps, DMRG_maxdim...)
    setcutoff!(sweeps, DMRG_cutoff...) 

    Parameters(N,t,U,ω,g0,g1,doping,sweeps,ceil(Int,N/2),T,τ,TEBD_cutoff,TEBD_maxdim)
end

function HubbardHolsteinModel(p::Parameters)
    N, t, U, ω, g0, g1, τ = p.N, p.t, p.U, p.ω, p.g0, p.g1, p.τ

    # make the sites 
    sites = siteinds("HubHolst", N) 

    # make the hamiltonian 
    ampo = Sum{Op}()
    for j=1:N-1
        # ∑_j,σ t * (c^†_jσ c_{j+1}σ + h.c.)
        ampo += -t,"Cdagup",j,"Cup",j+1
        ampo += -t,"Cdagup",j+1,"Cup",j
        ampo += -t,"Cdagdn",j,"Cdn",j+1
        ampo += -t,"Cdagdn",j+1,"Cdn",j
        
        # ∑_j U * n_j↑ n_j↓ 
        ampo += U,"Nupdn",j,"I",j

        # ∑_j ω * nb_j
        ampo += ω,"Nb",j

        # # ∑_j g0 * nf_j (b^†_j + b_j)
        ampo += g0,"Ntot(Bd+B)",j

        # # ∑_⟨ij⟩ g1 * nf_j (b^†_i + b_i)
        ampo += g1,"Ntot",j,"Bdag+B",j+1 
        ampo += g1,"Ntot",j+1,"Bdag+B",j
    end
    # Edge site
    ampo += U,"Nupdn",N
    ampo += ω,"Nb",N
    ampo += g0,"Ntot(Bd+B)",N
    
    H = MPO(ampo,sites)

    # make the trotter gates 
    gates = ITensor[]
    hgates = ITensor[]
    for j=1:N-1
        s1 = sites[j] # site j
        s2 = sites[j+1] # site j+1

        hj = -t*(op("Cdagup*F",s1) * op("Cup",s2)  # t * (c^†_jσ c_{j+1}σ + h.c.)
                 -op("Cup*F",s1) * op("Cdagup",s2) 
                 +op("Cdagdn*F",s1) * op("Cdn",s2) 
                 -op("Cdn*F",s1) * op("Cdagdn",s2)) 

            + U*(op("Nupdn",s1) * op("I",s2))    

            + ω*(op("Nb",s1) * op("I",s2))   

            + g0*(op("Ntot(Bd+B)",s1) * op("I",s2))

            + g1*(op("Ntot",s1) * op("Bdag+B",s2))
            + g1*(op("Bdag+B",s1) * op("Ntot",s2))
            

        Gj = exp(-1.0im * τ/2 * hj)
        push!(gates,Gj)
        push!(hgates, hj)
    end

    # End site 
    hn = U*op("Nupdn",sites[N]) + ω*op("Nb",sites[N]) + 
            g0*op("Ntot(Bd+B)",sites[N])
    Gn = exp(-1.0im * τ/2 * hn)
    push!(gates,Gn)
    push!(hgates, hn)

    append!(gates,reverse(gates))
    HubbardHolsteinModel(sites, H, gates, ampo)
end

function initialize_wavefcn(HH::HubbardHolsteinModel, p::Parameters)
    state = [isodd(n) ? "Up,0" : "Dn,0" for n=1:p.N] # NOTE: the QN of this state is preserved through DMRG

    # Account for doping
    if p.doping > 0
        spacing = floor(1/p.doping)
        state[1:spacing:end] .= "Emp,0"
    end

    productMPS(HH.sites,state) 
end

function run_DMRG(HH::HubbardHolsteinModel, p::Parameters)
    ϕ0 = initialize_wavefcn(HH,p)
    @show flux(ϕ0)
    energy, ϕ = dmrg(HH.mpo, ϕ0, p.sweeps)
    entropy = compute_entropy(ϕ, p.mid)
    return DMRGResults(ϕ, energy, entropy)
end

function apply_onesite_operator(ϕ::MPS, opname::String, sites, siteidx::Int)
    ψ = copy(ϕ) # Make a copy of the original state
    new_ψj = op(opname,sites[siteidx]) * ψ[siteidx] # Apply the local operator
    noprime!(new_ψj) 
    ψ[siteidx] = new_ψj
    return ψ
end

function compute_entropy(ψ::MPS, s)
    orthogonalize!(copy(ψ), s)
    _,S,_ = svd(ψ[s], (linkind(ψ, s-1), siteind(ψ,s)))
    S = S ./ sum(S.^2) # Normalize so that the sum squared of eigenvalues is 1
    SvN = 0.0
    for n=1:dim(S, 1)
        p = S[n,n]
        SvN -= p * log(p)
    end
    return SvN
end

function compute_overlap(ψ1::MPS, ψ2::MPS)
    LinearAlgebra.norm(inner(ψ1, ψ2))
end

function compute_phonon_number(ψ::MPS)
    expect(ψ,"Nb")
end

function time_evolve(dmrg_results::DMRGResults, HH::HubbardHolsteinModel, p::Parameters; order=4)
    ℋ = HH.ampo
    s = HH.sites
    ψ₀ = dmrg_results.ground_state
    nsteps = 10#floor(p.T/p.τ)
    t = 0.1
    println("Making 𝒰")
    𝒰 = exp(im * t * ℋ; alg=Trotter{order}(nsteps))
    println("Making U")
    U = Prod{ITensor}(𝒰, s)
    println("Applying the gate")
    ψ = apply(U,ψ₀)
    @show inner(ψ₀, ψ)

    #H = ITensor(ℋ, s)
    #𝒰ʳᵉᶠ = exp(im * t * ℋ)
    #Uʳᵉᶠ = ITensor(𝒰ʳᵉᶠ, s)
    #Uʳᵉᶠψ₀ = replaceprime(Uʳᵉᶠ * prod(ψ₀), 1 => 0)
end

function make_spectral_fcn(corrs, p::Parameters)
    # now has dimensions time x space
    sqw = fft.fftshift(fft.fft2(corrs * p.τ, norm="ortho"),axes=0) # TODO: Check the axis
    qs = range(0, stop=2*π, length=p.N+2)[2:end-1]
    ωs = 2 * π * fft.fftshift(fft.fftfreq(Int(p.T/p.τ), p.τ))
    for i in 1:p.N
        sqw[:,i] = imag.(exp(1im * qs[i] * p.mid) * sqw[:,i])
    end
    return real(sqw)/π, ωs, qs
end

function plot_correlation_function(tebd_results::TEBDResults)
    heatmap(LinearAlgebra.norm.(tebd_results.corrs'))
end

function plot_spectral_function(tebd_results::TEBDResults, p::Parameters; lims=nothing)
    ff, ωs, qs = make_spectral_fcn(tebd_results.corrs', p)
    if !isnothing(lims)
        ff = ff[lims[1]:lims[2],:]
        ωs = ωs[lims[1]:lims[2]]
    end
    ff = circshift(ff', p.mid-1)'
    qs = range(-π, stop=π, length=p.N+2)[2:end-1]
    maxval = maximum(abs.(ff))
    heatmap(qs, ωs, abs.(ff), c=:bwr, clims=(-maxval, maxval))
end

function plot_entropy(tebd_results::TEBDResults)
    niters = length(tebd_results.entropy)
    ϕ_entropy = [tebd_results.entropy[n][1] for n in 1:niters]
    ψ_entropy = [tebd_results.entropy[n][2] for n in 1:niters]
    plot(1:niters, ϕ_entropy, label="ϕ(t)")
    plot!(1:niters, ψ_entropy, label="ψ(t)")
    title!("Von Neumann Entropy")
    xlabel!("Iteration")
end

function plot_phonon_flux(tebd_results::TEBDResults)
    ## TO DO
end

function plot_overlap(tebd_results::TEBDResults)
    plot(1:length(tebd_results.self_overlap), tebd_results.self_overlap)
end

## CODE ##

## PARAMETERS ## 

# Model 
N = 8
t = 1 ## THIS TERM IS FINE
U = 8
ω = 0*t ## THIS TERM IS FINE (by itself)
g0 = 0*t ## THIS TERM IS FINE (by itself)
g1 = 0*g0 ## THIS TERM IS FINE (by itself)

# Simulation 
T = 5
τ = 0.01
DMRG_numsweeps = 20
TEBD_maxdim = 800
TEBD_cutoff = 1E-14

## TODO: Modify number of phonons on each site from script 
## TODO: Confirm that it's okay to overwrite the operators in the hubbardholstein.jl sites file
## TODO: How many phonons should we initially add to each site? 

# Specify operators of interest
A_t0 = "Cup"
A_t = "Cup"

# Initialize 
println("Initializing...")
params = parameters(N=N, t=t, U=U, ω=ω, g0=g0, g1=g1, 
                    DMRG_numsweeps=DMRG_numsweeps, 
                    T=T, τ=τ, TEBD_cutoff=TEBD_cutoff)
hubbholst = HubbardHolsteinModel(params)

# Run DMRG
println("Finding ground state...")
dmrg_results = run_DMRG(hubbholst, params)
@show compute_phonon_number(dmrg_results.ground_state)

# Time evolve
println("Time evolving...")
tebd_results = time_evolve(dmrg_results, hubbholst, params)

using Pkg; 
Pkg.activate("././.")
using Revise
using ITensors
using SciPy: fft
using Plots
using NPZ
using JLD

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
        DMRG_cutoff = 1E-10
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
    ampo = OpSum()
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
    for j=1:N-1
        s1 = sites[j] # site j
        s2 = sites[j+1] # site j+1

        hj_twosite = -t*(op("Cdagup*F",s1) * op("Cup",s2)  # t * (c^†_jσ c_{j+1}σ + h.c.)
                 -op("Cup*F",s1) * op("Cdagup",s2) 
                 +op("Cdagdn*F",s1) * op("Cdn",s2) 
                 -op("Cdn*F",s1) * op("Cdagdn",s2)) 
                + g1*(op("Ntot",s1) * op("Bdag+B",s2))
                + g1*(op("Bdag+B",s1) * op("Ntot",s2))

        hj_onesite = U*(op("Nupdn",s1) * op("I",s2))    
                    + ω*(op("Nb",s1) * op("I",s2))   
                    + g0*(op("Ntot(Bd+B)",s1) * op("I",s2))
            
        Gj_twosite = exp(-1.0im * τ/2 * hj_twosite)
        Gj_onesite = exp(-1.0im * τ/2 * hj_onesite)
        push!(gates,Gj_twosite)
        push!(gates,Gj_onesite)
    end
    # End site 
    hn = U*op("Nupdn",sites[N]) 
        + ω*op("Nb",sites[N]) 
        + g0*op("Ntot(Bd+B)",sites[N])
    Gn = exp(-1.0im * τ/2 * hn)
    push!(gates,Gn)
    append!(gates,reverse(gates))

    HubbardHolsteinModel(sites, H, gates)
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

function run_DMRG(HH::HubbardHolsteinModel, p::Parameters; alg="divide_and_conquer")
    ϕ0 = initialize_wavefcn(HH,p)
    @show flux(ϕ0)
    energy, ϕ = dmrg(HH.mpo, ϕ0, p.sweeps, alg=alg)
    entropy = compute_entropy(ϕ, p.mid)
    return DMRGResults(ϕ, energy, entropy)
end

function apply_twosite_operator(ϕ::MPS, opname_1::String, opname_2::String, sites, siteidx::Int)
    ϕ = copy(ϕ)

end

""" 
    function apply_onesite_operator(ϕ::MPS, opname::String, sites, siteidx::Int)
This is a wrapper function for the apply_op function, which actually applies the operator
Here, we just account for fermionic strings 
"""
function apply_onesite_operator(ϕ::MPS, opname::String, sites, siteidx::Int)
    ϕ = copy(ϕ) 

    ## Account for fermion sign using Jordan-Wigner strings ##
    if opname == "Cup" || opname == "Cdn"
        ϕ = apply_op(ϕ, op(opname,sites[siteidx]), siteidx)
        for i in reverse(1:(siteidx-1)) # Don't act with string on-site
            ϕ = apply_op(ϕ, op("F",sites[i]), i)
        end
        return ϕ
    elseif opname == "Cdagup" || opname == "Cdagdn"
        for i in 1:(siteidx-1) # Don't act with string on-site
            ϕ = apply_op(ϕ, op("F",sites[i]), i)
        end
        ϕ = apply_op(ϕ, op(opname,sites[siteidx]), siteidx)
        return ϕ
    end
    println("No Jordan-Wigner strings")

    # Otherwise, just apply the operator as usual
    return apply_op(ϕ, op(opname,sites[siteidx]), siteidx)
end

function apply_op(ϕ::MPS, op::ITensor, siteidx::Int)
    ϕ = copy(ϕ) # Make a copy of the original state
    orthogonalize!(ϕ, siteidx)
    new_ϕj = op * ϕ[siteidx] # Apply the local operator
    noprime!(new_ϕj) 
    ϕ[siteidx] = new_ϕj
    return ϕ
end

function compute_entropy(ψ::MPS, b::Int)
    orthogonalize!(ψ, b)
    U,S,V = svd(ψ[b], (linkind(ψ, b-1), siteind(ψ,b)))
    SvN = 0.0
    for n=1:dim(S, 1)
        p = S[n,n]^2
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

function equilibrium_correlations(dmrg_results::DMRGResults, corrtype::String, HH::HubbardHolsteinModel, p::Parameters)
    ϕ = copy(dmrg_results.ground_state)
    j = p.mid
    sites = HH.sites
    if corrtype=="spin"
        return correlation_matrix(ϕ, "Sz", "Sz")[:,j]
    elseif corrtype=="charge"
        ninj = correlation_matrix(ϕ, "Ntot", "Ntot")[:,j]
        ni = expect(ϕ, "Ntot")
        nj = ni[j]
        return ninj - nj .* ni
    elseif corrtype=="sSC"
        Σ_j = op("Cupdn", sites[j])
        ψ = apply_op(ϕ, Σ_j, j)
        function compute_corr(i::Int)
            Σ_i = op("Cdagupdn", sites[i])
            Σ_iψ = apply_op(ψ, Σ_i, i)
            return inner(ϕ,Σ_iψ)
        end
        return compute_corr.(collect(1:p.N))
    elseif corrtype=="pSC"
        Π_j = 1/sqrt(2)*(op("Cup",sites[j])*op("Cdn",sites[j+1]) + op("Cdn",sites[j])*op("Cup",sites[j+1]))
        ψ = apply_op(ϕ, Π_j, j)
        function compute_corr(i::Int)
            Π_i = 1/sqrt(2)*(op("Cdagdn",sites[i])*op("Cdagup",sites[i+1]) + op("Cdagup",sites[i])*op("Cdagdn",sites[i+1]))
            Π_iψ = apply_op(ψ, Π_i, i)
            return inner(ϕ,Π_iψ)
        end
        return compute_corr.(collect(1:p.N))
    elseif corrtype=="dSC"
        Δ_j = 1/sqrt(2)*(op("Cup",sites[j])*op("Cdn",sites[j+1]) - op("Cdn",sites[j])*op("Cup",sites[j+1]))
        ψ = apply_op(ϕ, Δ_j, j)
        function compute_corr(i::Int)
            Δ_i = 1/sqrt(2)*(op("Cdagdn",sites[i])*op("Cdagup",sites[i+1]) - op("Cdagup",sites[i])*op("Cdagdn",sites[i+1]))
            Δ_iψ = apply_op(ψ, Δ_i, i)
            return inner(ϕ,Δ_iψ)
        end
        return compute_corr.(collect(1:p.N))
    else
        @error "Unknown spin correlation type"
    end
end

function plot_equilibrium_correlations(dmrg_results::DMRGResults, corrtype::String, HH::HubbardHolsteinModel, p::Parameters)
    corrs = equilibrium_correlations(dmrg_results,corrtype,HH,p)
    plot(1:length(corrs), corrs)
    title!(corrtype)
    xlabel!("Site")
    ylabel!("Correlation")
end

function compute_correlations(dmrg_results::DMRGResults, A_t0::String, A_t::String, HH::HubbardHolsteinModel, p::Parameters)
    # Results 
    corrs = []
    entropy = []
    self_overlap = []
    phonon_flux = []

    # The wavefunction being acted upon at t=0, |ψ⟩ = A_t0|ϕ⟩
    ϕ = copy(dmrg_results.ground_state)
    ψ = copy(ϕ)
    
    # Apply A_t0 to middle site
    ψ = apply_onesite_operator(ψ, A_t0, HH.sites, p.mid)

    # Measure correlation before time evolution
    A_tψ = apply_onesite_operator(ψ, A_t, HH.sites, p.mid)
    # Measure the correlation fcn 
    @show inner(ϕ, A_tψ)

    # Parameters for time evolution
    nsteps = floor(p.T/p.τ) # Number of time steps for time evolution
    t = 0.0
    for step in 1:nsteps
        print(floor(Int,step),"-")
        ϕ = apply(HH.gates, ϕ; maxdim=p.TEBD_maxdim, cutoff=p.TEBD_cutoff)
        ψ = apply(HH.gates, ψ; maxdim=p.TEBD_maxdim, cutoff=p.TEBD_cutoff) # evolve forward

        t += p.τ 

        ### SANITY CHECKS
        if step%1==0
            # Compute entropy
            push!(entropy, (compute_entropy(ϕ, p.mid),compute_entropy(ψ, p.mid)))
            println("Entropy of ϕ, ψ: ", entropy[end])

            # Compute overlap of ϕ with its original self
            push!(self_overlap, compute_overlap(ϕ,dmrg_results.ground_state))
            println("Overlap of ϕ with itself: ", self_overlap[end])

            # Phonon flux
            push!(phonon_flux, (compute_phonon_number(ϕ), compute_phonon_number(ψ)))
            println("Phonon flux ϕ: ", sum(phonon_flux[end][1])/p.N)
            println("Phonon flux ψ: ", sum(phonon_flux[end][2])/p.N)
        end

        # Calculate ⟨ϕ(t)|c_j^† c_i|ϕ(0)⟩
        function measure_corr(j::Int)
            # Apply the second measurement operator
            A_tψ = apply_onesite_operator(ψ, A_t, HH.sites, j)
            return dot(ϕ,A_tψ)
        end

        # Measure the correlation fcn 
        push!(corrs,measure_corr.(collect(1:p.N)))
    end
    return TEBDResults(entropy, self_overlap, phonon_flux, hcat(corrs...))
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

function compare_to_ED(tebd_results::TEBDResults, p::Parameters)
    @assert p.τ == 0.01 # Time scale of the ED results
    corrs = tebd_results.corrs
    ED_corrs = npzread("/Users/nicole/Dropbox/Grad school/Devereaux lab/ITensors.jl/examples/dmrg/ed.npy")
    numsteps = size(corrs)[2]
    plot(1:numsteps, real.(corrs[4,:]), label="DMRG real part")
    plot!(1:numsteps, real.(ED_corrs[1:numsteps]), label="ED real part")
    plot!(1:numsteps, imag.(corrs[4,:]), label="DMRG complex part")
    plot!(1:numsteps, imag.(ED_corrs[1:numsteps]), label="ED complex part")
end

function plot_phonon_flux(tebd_results::TEBDResults)
    ## TO DO
end

function plot_overlap(tebd_results::TEBDResults)
    plot(1:length(tebd_results.self_overlap), tebd_results.self_overlap)
end

function save_structs(struc, path::String)
    function Name(arg)
        string(arg)
    end
    fnames = fieldnames(typeof(struc))
    for fn in fnames 
        n = Name(fn)
        d = getfield(struc, fn)
        jldopen(path, "w") do file
            write(file, n, d) 
        end
    end
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
T = 10#40
τ = 0.01#0.05
DMRG_numsweeps = 40
DMRG_maxdim = 600
TEBD_maxdim = 800
TEBD_cutoff = 1E-14
DMRG_cutoff = 1E-14

# Saveout info 
save_path = "/Users/nicole/Dropbox/Grad school/Devereaux lab/ITensors.jl/examples/dmrg/results/1dhubbholst_output.jld"

## TODO: Modify number of phonons on each site from script 

# Specify operators of interest
A_t0 = "Cup"
A_t = "Cdagup"

# Initialize 
println("Initializing...")
params = parameters(N=N, t=t, U=U, ω=ω, g0=g0, g1=g1, 
                    DMRG_numsweeps=DMRG_numsweeps,
                    DMRG_maxdim=DMRG_maxdim, DMRG_cutoff=DMRG_cutoff,
                    T=T, τ=τ, TEBD_cutoff=TEBD_cutoff)
hubbholst = HubbardHolsteinModel(params)
save_structs(params, save_path)

# Run DMRG
println("Finding ground state...")
dmrg_results = run_DMRG(hubbholst, params, alg="divide_and_conquer")
save_structs(dmrg_results, save_path)

# Compute correlation functions 
println("Computing correlation functions...")
tebd_results = compute_correlations(dmrg_results, A_t0, A_t, hubbholst, params)
save_structs(tebd_results, save_path)
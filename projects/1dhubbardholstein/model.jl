using ITensors
include("dmrg_lbo.jl")

## STRUCTS ##

struct Parameters
    # Model parameters
    N::Int
    t::Real
    U::Real
    ω::Real
    g0::Real
    g1::Real
    λ::Real
    doping::Real
    max_phonons::Int
    init_phonons::Int

    # DMRG parameters
    DMRG_numsweeps
    DMRG_noise
    DMRG_maxdim
    DMRG_cutoff
    DMRG_LBO::Bool
    max_LBO_dim::Int
    min_LBO_dim::Int

    # TEBD parameters
    mid::Int
    T::Int
    τ::Real
    TEBD_cutoff
    TEBD_maxdim
    TEBD_LBO::Bool
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
    ground_state_entropy 
    optimized_basis
end

struct EquilibriumCorrelations
    spin
    charge
    sSC
    pSC
    dSC
end

struct TEBDResults
    entropy
    self_overlap
    corrs 
    phi_t 
    psi_t
end

## SETTING UP THE MODEL ## 

function parameters(;N::Int, t::Real, U::Real=nothing, ω::Real=nothing, 
                    g0::Real=nothing, g1::Real=nothing, λ::Real=nothing, 
                    doping::Real=0, max_phonons::Int=1, init_phonons::Int=0,
                    DMRG_numsweeps::Int=20, DMRG_noise=nothing, 
                    DMRG_maxdim=nothing, DMRG_cutoff=nothing, DMRG_LBO=false,
                    max_LBO_dim=nothing, min_LBO_dim=4,
                    T::Int=25, τ::Real=0.1, TEBD_cutoff=1E-14, TEBD_maxdim=400, TEBD_LBO=false, )
    
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
        g1 = 0
    end
    if isnothing(λ)
        λ = 0
    end
    if max_phonons==0
        @assert ω==g0==g1==0
    end
    if ω==0
        @assert λ==0
    end

    if isnothing(DMRG_noise)
        DMRG_noise = [1E-6,1E-6,1E-8,0]
    end
    if isnothing(DMRG_maxdim)
        DMRG_maxdim = [20,40,100,200,400]
    end
    if isnothing(DMRG_cutoff)
        DMRG_cutoff = 1E-10
    end
    if isnothing(max_LBO_dim)
        max_LBO_dim = 12
    end

    Parameters(N,t,U,ω,g0,g1,λ,doping,max_phonons,init_phonons,
                DMRG_numsweeps,DMRG_noise,DMRG_maxdim,DMRG_cutoff,DMRG_LBO,
                max_LBO_dim,min_LBO_dim,
                ceil(Int,N/2),T,τ,TEBD_cutoff,TEBD_maxdim,TEBD_LBO)
end

"""
Use this function for reconstructing from scratch given site indices 
"""
function HubbardHolsteinModel(p::Parameters, sites::Vector{Index{Vector{Pair{QN, Int64}}}})
    N, t, U, ω, g0, g1, λ, τ, max_phonons = p.N, p.t, p.U, p.ω, p.g0, p.g1, p.λ, p.τ, p.max_phonons

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
        if max_phonons >= 1
            ampo += ω,"Nb",j

            # # ∑_j g0 * nf_j (b^†_j + b_j)
            ampo += g0,"Ntot(Bd+B)",j

            # # ∑_⟨ij⟩ g1 * nf_j (b^†_i + b_i)
            ampo += g1,"Ntot",j,"Bdag+B",j+1
            ampo += g1,"Ntot",j+1,"Bdag+B",j

            # quartic term
            ampo += λ,"Nb^2",j
        end
    end
    # Edge site
    ampo += U,"Nupdn",N
    if max_phonons >= 1
        ampo += ω,"Nb",N
        ampo += g0,"Ntot(Bd+B)",N
        ampo += λ,"Nb^2",N
    end
    H = MPO(ampo,sites)

    ## Helper functions for making the Trotter gates ##
    function make_twosite(s1, s2)
        hj_twosite = -t*(op("Cdagup*F",s1) * op("Cup",s2)  # t * (c^†_jσ c_{j+1}σ + h.c.)
                 -op("Cup*F",s1) * op("Cdagup",s2) 
                 +op("Cdagdn*F",s1) * op("Cdn",s2) 
                 -op("Cdn*F",s1) * op("Cdagdn",s2)) 
        if max_phonons >= 1
            hj_twosite  = hj_twosite
                + g1*(op("Ntot",s1) * op("Bdag+B",s2))
                + g1*(op("Bdag+B",s1) * op("Ntot",s2))
        end

        return exp(-1.0im * τ/2 * hj_twosite)
    end

    function make_onesite(s1)
        hj_onesite = U * op("Nupdn",s1)
        if max_phonons >=1
            hj_onesite = hj_onesite
                    + ω * op("Nb",s1)   
                    + g0 * op("Ntot(Bd+B)",s1)
                    + λ * op("Nb^2",s1)
        end

        return exp(-1.0im * τ/2 * hj_onesite)
    end

    function make_endsite(sn)
        hn = U*op("Nupdn",sn) 
        if max_phonons >= 1
            hn = hn
            + ω*op("Nb",sn) 
            + g0*op("Ntot(Bd+B)",sn)
            + λ*op("Nb^2",sn) 
        end
        return exp(-1.0im * τ/2 * hn)
    end

    # make the trotter gates 
    gates = ITensor[]
    for j=1:N-1
        s1 = sites[j] # site j
        s2 = sites[j+1] # site j+1

        hj_twosite = make_twosite(s1, s2)
        hj_onesite = make_onesite(s1)
            
        push!(gates,hj_twosite)
        push!(gates,hj_onesite)
    end
    # End site 
    hn = make_endsite(sites[N])
    push!(gates,hn)

    # Add all the gates in reverse 
    append!(gates,reverse(gates))

    # Return the struct 
    HubbardHolsteinModel(sites, H, gates)
end

function HubbardHolsteinModel(p::Parameters)
    # make the sites 
    sites = siteinds("HubHolst", p.N; dim=p.max_phonons+1)

    return HubbardHolsteinModel(p, sites)
end

## GENERIC STATE FUNCTIONS ##

function compute_overlap(ψ1::MPS, ψ2::MPS)
    LinearAlgebra.norm(inner(ψ1, ψ2))
end

function compute_phonon_number(ψ::MPS)
    expect(ψ,"Nb")
end

function compute_electron_number(ψ::MPS)
    expect(ψ,"Ntot")
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

## DMRG ## 

function initialize_wavefcn(HH::HubbardHolsteinModel, p::Parameters)
    # Extract the information 
    ups = "Up,"*string(p.init_phonons)
    downs = "Dn,"*string(p.init_phonons)
    emps = "Emp,"*string(p.init_phonons)
    
    state = [isodd(n) ? ups : downs for n=1:p.N] 

    # Account for doping
    if p.doping > 0
        spacing = floor(Int,1/p.doping)
        state[1:spacing:end] .= emps
    end

    # NOTE: the QN of this state is preserved during DMRG
    productMPS(HH.sites,state) 
end

function run_DMRG(HH::HubbardHolsteinModel, p::Parameters; alg="divide_and_conquer")
    # Set DMRG params
    sweeps = Sweeps(p.DMRG_numsweeps)
    setnoise!(sweeps, p.DMRG_noise...) # Very important to use noise for this model
    setmaxdim!(sweeps, p.DMRG_maxdim...)
    setcutoff!(sweeps, p.DMRG_cutoff...) 
    
    ϕ0 = initialize_wavefcn(HH,p)
    @show flux(ϕ0)
    if p.DMRG_LBO # If performing local basis optimization
        energy, ϕ, Rs = dmrg_lbo(HH.mpo, ϕ0, sweeps, alg=alg, LBO=true, 
                                    max_LBO_dim=p.max_LBO_dim, min_LBO_dim=p.min_LBO_dim)
    else
        energy, ϕ = dmrg(HH.mpo, ϕ0, sweeps, alg=alg)
        Rs = nothing
    end
    entropy = compute_entropy(ϕ, p.mid)
    return DMRGResults(ϕ, energy, entropy, Rs)
end

## CORRELATION FUNCTIONS ## 

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

    # Otherwise, just apply the operator as usual
    return apply_op(ϕ, op(opname,sites[siteidx]), siteidx)
end

function apply_twosite_operator(ϕ::MPS, opname::String, sites, siteidx::Int)
    ϕ = copy(ϕ)
    if opname == "pSC"
        # Apply F string first 
        ϕ = apply_op(ϕ, op("F",sites[siteidx]), siteidx)

        # Then apply the operator
        Π = 1/sqrt(2)*(op("Cup",sites[siteidx])*op("Cdn",sites[siteidx+1])
                            + op("Cdn",sites[siteidx])*op("Cup",sites[siteidx+1]))
        return apply_op_twosite(ϕ,Π,siteidx)
                    
    elseif opname == "pSCdag"
        # Apply the operator first
        Π_dag = 1/sqrt(2)*(op("Cdagdn",sites[siteidx+1])*op("Cdagup",sites[siteidx])
                            + op("Cdagup",sites[siteidx+1])*op("Cdagdn",sites[siteidx]))
        ϕ = apply_op_twosite(ϕ,Π_dag,siteidx)

        # Apply the F string 
        return apply_op(ϕ, op("F",sites[siteidx]), siteidx)

    elseif opname == "dSC"
        # Apply F string first 
        ϕ = apply_op(ϕ, op("F",sites[siteidx]), siteidx)

        # Then apply the operator
        Δ = 1/sqrt(2)*(op("Cup",sites[siteidx])*op("Cdn",sites[siteidx+1])
                            - op("Cdn",sites[siteidx])*op("Cup",sites[siteidx+1]))
        return apply_op_twosite(ϕ,Δ,siteidx)

    elseif opname == "dSCdag"
        # Apply the operator first
        Δ_dag = 1/sqrt(2)*(op("Cdagdn",sites[siteidx+1])*op("Cdagup",sites[siteidx])
                            - op("Cdagup",sites[siteidx+1])*op("Cdagdn",sites[siteidx]))
        ϕ = apply_op_twosite(ϕ,Δ_dag,siteidx)

        # Apply the F string 
        return apply_op(ϕ, op("F",sites[siteidx]), siteidx)
    end
    @error "No recognized two-site operator"
end

function apply_op(ϕ::MPS, op::ITensor, siteidx::Int)
    ϕ = copy(ϕ) # Make a copy of the original state
    orthogonalize!(ϕ, siteidx)
    new_ϕj = op * ϕ[siteidx] # Apply the local operator
    noprime!(new_ϕj) 
    ϕ[siteidx] = new_ϕj
    return ϕ
end

function apply_op_twosite(ϕ::MPS, G::ITensor, siteidx::Int; cutoff=1E-8)
    ϕ = copy(ϕ)
    # Note: siteidx corresponds to the leftermost site
    orthogonalize!(ϕ,siteidx)
    wf = (ϕ[siteidx] * ϕ[siteidx+1]) * G
    noprime!(wf)

    inds_site = uniqueinds(ϕ[siteidx],ϕ[siteidx+1])
    U,S,V = svd(wf,inds_site,cutoff=cutoff)
    ϕ[siteidx] = U
    ϕ[siteidx+1] = S*V
    return ϕ
end

function compute_all_equilibrium_correlations(dmrg_results::DMRGResults, 
                                            HH::HubbardHolsteinModel;
                                            start=nothing, stop=nothing)

    N = length(HH.sites)
    if isnothing(start)
        start = floor(Int,0.25*N)
    end 
    if isnothing(stop)
        stop = ceil(Int,0.75*N)
    end

    corrtypes = ["spin","charge","sSC","pSC","dSC"]
    corrs = []
    for corrtype in corrtypes
        corr = equilibrium_correlations(dmrg_results,corrtype,HH,start,stop)
        push!(corrs,corr)
    end
    EquilibriumCorrelations(corrs...)
end

function equilibrium_correlations(dmrg_results::DMRGResults, corrtype::String, 
                                HH::HubbardHolsteinModel, start::Int, stop::Int,)
    
    ϕ = copy(dmrg_results.ground_state)
    j = start
    sites = HH.sites
    if corrtype=="spin"
        return correlation_matrix(ϕ, "Sz", "Sz")[start:stop,j]
    elseif corrtype=="charge"
        ninj = correlation_matrix(ϕ, "Ntot", "Ntot")[start:stop,j]
        ni = expect(ϕ, "Ntot")
        nj = ni[j]
        return ninj - nj .* (ni[start:stop])
    elseif corrtype=="sSC"
        ψ = apply_onesite_operator(ϕ, "Cupdn", sites, j)
        function compute_corr_sSC(i::Int)
            Σ_iψ = apply_onesite_operator(ψ, "Cdagupdn", sites, i)
            return inner(ϕ,Σ_iψ)
        end
        return compute_corr_sSC.(collect(start:stop))
    elseif corrtype=="pSC"
        ψ = apply_twosite_operator(ϕ, "pSC", sites, j)
        function compute_corr_pSC(i::Int)
            Π_iψ = apply_twosite_operator(ψ, "pSCdag", sites, i)
            return inner(ϕ,Π_iψ)
        end
        return compute_corr_pSC.(collect(start:stop))
    elseif corrtype=="dSC"
        ψ = apply_twosite_operator(ϕ, "dSC", sites, j)
        function compute_corr_dSC(i::Int)
            Δ_iψ = apply_twosite_operator(ψ, "dSCdag", sites, i)
            return inner(ϕ,Δ_iψ)
        end
        return compute_corr_dSC.(collect(start:stop))        
    end
    @error "Unknown spin correlation type"
end

function compute_correlations(dmrg_results::DMRGResults, 
    A_t0::String, A_t::String, 
    HH::HubbardHolsteinModel, p::Parameters;
    interim_save=false, savepath=nothing)

    if p.TEBD_LBO
        return tebd_lbo_corrs(dmrg_results, 
        A_t0, A_t, 
        HH, p;
        interim_save=interim_save, savepath=savepath)
    else
        return tebd_corrs(dmrg_results, 
        A_t0, A_t, 
        HH, p;
        interim_save=interim_save, savepath=savepath)
    end
end

function tebd_lbo_corrs(dmrg_results::DMRGResults, 
                        A_t0::String, A_t::String, 
                        HH::HubbardHolsteinModel, p::Parameters;
                        interim_save=false, savepath=nothing)

    # The wavefunction being acted upon at t=0, |ψ⟩ = A_t0|ϕ⟩
    ϕ = copy(dmrg_results.ground_state)
    ψ = copy(ϕ)
    Rs = copy(dmrg_results.optimized_basis)
    gates = copy(HH.gates)
    
    # Apply A_t0 to middle site
    ψ = apply_onesite_operator(ψ, A_t0, HH.sites, p.mid)

    # Measure correlation before time evolution
    A_tψ = apply_onesite_operator(ψ, A_t, HH.sites, p.mid)
    # Measure the correlation fcn 
    @show inner(ϕ, A_tψ)

    # Parameters for time evolution
    nsteps = floor(p.T/p.τ) # Number of time steps for time evolution
    t = 0.0

    # Results 
    corrs = []
    entropy = []
    self_overlap = Float64[]

    for step in 1:nsteps
        print(floor(Int,step),"-")

        ## TO DO: IMPLEMENT LBO FOR THE TEBD STEP ## 
        # 1. Take the gate and contract (on only ONE side) with the optimized basis   
        # 2. Put ψ in the optimized basis 
        # 3. Contract ψ and gates in the optimized basis. Will be in the bare basis on the uncontracted side
        # 4. Optimize the basis of ψ again using the local reduced density matrix 
        # 5. Update the rotation matrix 

        ϕ = apply(HH.gates, ϕ; maxdim=p.TEBD_maxdim, cutoff=p.TEBD_cutoff)
        ψ = apply(HH.gates, ψ; maxdim=p.TEBD_maxdim, cutoff=p.TEBD_cutoff) 

        t += p.τ 

        ### SANITY CHECKS
        if step%1==0
            # Compute entropy
            push!(entropy, [compute_entropy(ϕ, p.mid),compute_entropy(ψ, p.mid)])
            println("Entropy of ϕ, ψ: ", entropy[end])

            # Compute overlap of ϕ with its original self
            push!(self_overlap, compute_overlap(ϕ,dmrg_results.ground_state))
            println("Overlap of ϕ with itself: ", self_overlap[end])
        end

        # Calculate ⟨ϕ(t)|c_j^† c_i|ϕ(0)⟩
        function measure_corr(j::Int)
            # Apply the second measurement operator
            A_tψ = apply_onesite_operator(ψ, A_t, HH.sites, j)
            return dot(ϕ,A_tψ)
        end

        # Measure the correlation fcn 
        push!(corrs,measure_corr.(collect(1:p.N)))

        if interim_save && step%10==0
            @assert !isnothing(savepath)
            tebd_results_interim = TEBDResults(hcat(entropy...), self_overlap, 
                                hcat(corrs...), ϕ, ψ)
            save_structs(tebd_results_interim, save_path)
        end
            
    end

    return TEBDResults(hcat(entropy...), self_overlap, 
                        hcat(corrs...), ϕ, ψ)

end

function tebd_corrs(dmrg_results::DMRGResults, 
                            A_t0::String, A_t::String, 
                            HH::HubbardHolsteinModel, p::Parameters;
                            interim_save=false, savepath=nothing)

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

    # Results 
    corrs = []
    entropy = []
    self_overlap = Float64[]

    for step in 1:nsteps
        print(floor(Int,step),"-")

        ϕ = apply(HH.gates, ϕ; maxdim=p.TEBD_maxdim, cutoff=p.TEBD_cutoff)
        ψ = apply(HH.gates, ψ; maxdim=p.TEBD_maxdim, cutoff=p.TEBD_cutoff) 

        t += p.τ 

        ### SANITY CHECKS
        if step%1==0
            # Compute entropy
            push!(entropy, [compute_entropy(ϕ, p.mid),compute_entropy(ψ, p.mid)])
            println("Entropy of ϕ, ψ: ", entropy[end])

            # Compute overlap of ϕ with its original self
            push!(self_overlap, compute_overlap(ϕ,dmrg_results.ground_state))
            println("Overlap of ϕ with itself: ", self_overlap[end])
        end

        # Calculate ⟨ϕ(t)|c_j^† c_i|ϕ(0)⟩
        function measure_corr(j::Int)
            # Apply the second measurement operator
            A_tψ = apply_onesite_operator(ψ, A_t, HH.sites, j)
            return dot(ϕ,A_tψ)
        end

        # Measure the correlation fcn 
        push!(corrs,measure_corr.(collect(1:p.N)))

        if interim_save && step%10==0
            @assert !isnothing(savepath)
            tebd_results_interim = TEBDResults(hcat(entropy...), self_overlap, 
                                hcat(corrs...), ϕ, ψ)
            save_structs(tebd_results_interim, save_path)
        end
            
    end

    return TEBDResults(hcat(entropy...), self_overlap, 
                        hcat(corrs...), ϕ, ψ)
end



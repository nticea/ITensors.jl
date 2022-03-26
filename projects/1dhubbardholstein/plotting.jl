## IMPORTS ## 
using SciPy: fft
using Plots
using DSP
include(joinpath(@__DIR__,"model.jl"))

## EQUILIBRIUM CORRELATIONS ##

function plot_density(dmrg_results::DMRGResults)
    n = expect(dmrg_results.ground_state, "Ntot")
    plot(1:length(n), n)
    ylims!((0.2,1.2))
    ylabel!("⟨n⟩")
    xlabel!("Site")
    title!("Electron density")
end

function plot_phonon_flux(dmrg_results::DMRGResults)
    n = expect(dmrg_results.ground_state, "Nb")
    plot(1:length(n), n)
    ylims!((0.2,1.2))
    ylabel!("⟨nb⟩")
    xlabel!("Site")
    title!("Phonon density")
end

function plot_equilibrium_correlations(dmrg_results::DMRGResults, corrtype::String, 
                                        HH::HubbardHolsteinModel;
                                        start=nothing, stop=nothing)

    N = length(HH.sites)
    if isnothing(start)
        start = floor(Int, N/4)
    end
    if isnothing(stop)
        stop = ceil(Int, N/4*3)
    end

    corrs = equilibrium_correlations(dmrg_results,corrtype,HH,start,stop)

    plot(log10.(collect(0:stop-start)), log10.(abs.(corrs)))
    title!(corrtype*"-"*corrtype*" correlation")
    xlabel!("Distance from centre site (log10)")
    ylabel!("Correlation (log10)")
end

## CHECKING TEBD RESULTS ##

function plot_entropy(tebd_results::TEBDResults)
    ent = tebd_results.entropy
    niters = size(ent)[2]
    ϕ_entropy = ent[1,:]
    ψ_entropy = ent[2,:]
    plot(1:niters, ϕ_entropy, label="ϕ(t)")
    plot!(1:niters, ψ_entropy, label="ψ(t)")
    title!("Von Neumann Entropy")
    xlabel!("Iteration")
end

function plot_overlap(tebd_results::TEBDResults)
    plot(1:length(tebd_results.self_overlap), tebd_results.self_overlap)
end

## TIME-DEPENDENT CORRELATIONS ## 

function plot_correlation_function(tebd_results::TEBDResults)
    heatmap(LinearAlgebra.norm.(tebd_results.corrs'),color_palette=:dense)
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

function convolve(u;shape="exponential")
    if shape=="exponential"
        t = collect(1:(length(u)*2-1))
        λ = 1/sqrt(length(t))
        v = exp.(-λ*t)
    else
        @error "Not implemented"
    end
    u = convert(Array{ComplexF64,1}, u)
    return conv(u,v)
end

function plot_spectral_function(tebd_results::TEBDResults, p::Parameters; 
                                lims=nothing, do_convolve=false)
    # Optionally convolve raw time data with decaying exponential
    corrs = tebd_results.corrs' # num_time_steps x num_sites
    if do_convolve
        @error "Not implemented yet, sorry"
        # for (i,r) in enumerate(eachrow(corrs_raw))
        #     @show size(convolve(r))
        #     corrs[:,i] .= convolve(r)
        # end
    end
    ff, ωs, qs = make_spectral_fcn(corrs, p)
    if !isnothing(lims)
        ff = ff[lims[1]:lims[2],:]
        ωs = ωs[lims[1]:lims[2]]
    end
    ff = circshift(ff', p.mid-1)'
    qs = range(-π, stop=π, length=p.N+2)[2:end-1]
    maxval = maximum(abs.(ff))
    heatmap(qs, ωs, abs.(ff), yflip=true, c=:bwr, clims=(-maxval, maxval))
end

function plot_spectral_function_slice(tebd_results::TEBDResults, p::Parameters; slice=0)
    findnearest(A::AbstractArray,t) = findmin(abs.(A.-t))[2]
    
    ff, ωs, qs = make_spectral_fcn(tebd_results.corrs', p)
    ff = circshift(ff', p.mid-1)'
    qs = range(-π, stop=π, length=p.N+2)[2:end-1]
    ω = findnearest(ωs,slice)
    plot(qs, abs.(ff)[ω,:])
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
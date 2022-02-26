using Pkg; 
Pkg.activate("././.")
using Revise
using ITensors

#
# DMRG calculation of the extended Hubbard model
# ground state wavefunction, and spin densities
#

let
  ## PARAMS ##

  # Model
  N = 8
  t = 1.0
  U = 8

  # DMRG 
  sweeps = Sweeps(20)
  setmaxdim!(sweeps, 50, 100, 200, 400, 800, 800)
  setcutoff!(sweeps, 1E-12)
  @show sweeps

  # TEBD
  T = 10
  τ = 0.01
  TEBD_maxdim = 800
  TEBD_cutoff = 1E-14

  ## SET UP THE MODEL ## 

  # Make N electron sites with conserved number and spin 
  sites = siteinds("Electron", N; conserve_qns=true)

  # Make the Hamiltonian 
  ampo = OpSum()
  for i in 1:(N - 1)
    ampo += -t, "Cdagup", i, "Cup", i + 1
    ampo += -t, "Cdagup", i + 1, "Cup", i
    ampo += -t, "Cdagdn", i, "Cdn", i + 1
    ampo += -t, "Cdagdn", i + 1, "Cdn", i

    ampo += U, "Nupdn", i
  end
  # Edge site
  ampo += U, "Nupdn", N
  H = MPO(ampo, sites)

  ## FIND GROUND STATE ##
  
  # Initialize wavefunction 
  psi0 = productMPS(sites, n -> isodd(n) ? "Up" : "Dn")
  @show flux(psi0)

  # Start DMRG calculation
  energy, psi = dmrg(H, psi0, sweeps)

  ## TIME-EVOLVE ##

  # Make the trotter gates 
  gates = ITensor[]
  for j=1:N-1
    s1 = sites[j] # site j
    s2 = sites[j+1] # site j+1

    hj = -t * (op("Cdagup*F",s1) * op("Cup",s2)  # t * (c^†_jσ c_{j+1}σ + h.c.)
                 -op("Cup*F",s1) * op("Cdagup",s2) 
                 +op("Cdagdn*F",s1) * op("Cdn",s2) 
                 -op("Cdn*F",s1) * op("Cdagdn",s2)) 

    Gj = exp(-1.0im * τ/2 * hj)
    push!(gates,Gj)
    hj2 = U * op("Nupdn",s1)   
    Gj2 = exp(-1.0im * τ/2 * hj2)
    push!(gates,Gj2)
  end
  # End site 
  hn = U*op("Nupdn",sites[N]) 
  Gn = exp(-1.0im * τ/2 * hn)
  push!(gates,Gn)
  # Append gates in reverse to complete Trotter formula
  append!(gates,reverse(gates))

  function compute_overlap(ψ1::MPS, ψ2::MPS)
      LinearAlgebra.norm(inner(ψ1, ψ2))
  end

  # Evolve
  nsteps = floor(T/τ) # Number of time steps for time evolution
  t = 0.0
  ϕ = copy(psi)
  for step in 1:nsteps
    print(floor(Int,step),"-")
    ϕ = apply(gates, ϕ; maxdim=TEBD_maxdim, cutoff=TEBD_cutoff)
    t += τ 

    @show compute_overlap(ϕ, psi)

  end

end
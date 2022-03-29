using LinearAlgebra
using SparseArrays

"""
    space(::SiteType"HubHolst"; 
          conserve_qns = false,
          conserve_sz = conserve_qns,
          conserve_nf = conserve_qns,
          conserve_nfparity = conserve_qns,
          qnname_sz = "Sz",
          qnname_nf = "Nf",
          qnname_nfparity = "NfParity")

Create the Hilbert space for a site of type "HubHolst".

Optionally specify the conserved symmetries and their quantum number labels.
"""
function ITensors.space(
  ::SiteType"HubHolst";
  dim=2, #Default is 1 phonon
  conserve_qns=true,
  conserve_sz=conserve_qns,
  conserve_nf=conserve_qns,
  conserve_nfparity=conserve_qns,
  qnname_sz="Sz",
  qnname_nf="Nf",
  qnname_nfparity="NfParity")

  if conserve_sz && conserve_nf
    return [
      QN((qnname_nf, 0, -1), (qnname_sz, 0)) => 1 * dim 
      QN((qnname_nf, 1, -1), (qnname_sz, +1)) => 1 * dim
      QN((qnname_nf, 1, -1), (qnname_sz, -1)) => 1 * dim
      QN((qnname_nf, 2, -1), (qnname_sz, 0)) => 1 * dim
    ]
  elseif conserve_nf
    return [
      QN(qnname_nf, 0, -1) => 1 * dim
      QN(qnname_nf, 1, -1) => 2 * dim
      QN(qnname_nf, 2, -1) => 1 * dim
    ]
  elseif conserve_sz
    return [
      QN((qnname_sz, 0), (qnname_nfparity, 0, -2)) => 1 * dim
      QN((qnname_sz, +1), (qnname_nfparity, 1, -2)) => 1 * dim
      QN((qnname_sz, -1), (qnname_nfparity, 1, -2)) => 1 * dim
      QN((qnname_sz, 0), (qnname_nfparity, 0, -2)) => 1 * dim
    ]
  elseif conserve_nfparity
    return [
      QN(qnname_nfparity, 0, -2) => 1 * dim
      QN(qnname_nfparity, 1, -2) => 2 * dim
      QN(qnname_nfparity, 0, -2) => 1 * dim
    ]
  end
  return 4
end

function ITensors.val(::ValName{N}, ::SiteType"HubHolst") where {N}
    hubbard_type, phonon_num = split(String(N), ",")
    if hubbard_type == "Emp"
        n1 = 1
    elseif hubbard_type == "Up"
        n1 = 2
    elseif hubbard_type == "Dn"
        n1 = 3
    elseif hubbard_type == "UpDn"
        n1 = 4
    else
        throw(DomainError(hubbard_type, "expects Emp, Up, Dn, UpDn"))
    end

    n2 = parse(Int, String(phonon_num)) + 1
    return n1+n2
end

function ITensors.state(n::StateName{N}, ::SiteType"HubHolst", s::Index) where {N}
    hubbard_type, n = split(String(name(n)), ",")
    n = parse(Int, n)

    if hubbard_type == "Emp"
        st_hubb = [1, 0, 0, 0]
    elseif hubbard_type == "Up"
        st_hubb = [0, 1, 0, 0]
    elseif hubbard_type == "Dn"
        st_hubb = [0, 0, 1, 0]
    elseif hubbard_type == "UpDn"
        st_hubb = [0, 0, 0, 1]
    else
        throw(DomainError(hubbard_type, "expects Emp, Up, Dn, UpDn"))
    end

    st_ph = zeros(floor(Int,dim(s)/4))
    st_ph[n+1] = 1

    return kron(st_hubb, st_ph)
end

## ELECTRON OPERATORS ## 

function ITensors.op(::OpName"I", ::SiteType"HubHolst", s::Index)
  M = Matrix(I, Int(dim(s)), Int(dim(s)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Nup", ::SiteType"HubHolst", s::Index)
  Nup = [
    0.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 1.0
  ]
  M = kron(Nup, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Ndn", ::SiteType"HubHolst", s::Index)
  Ndn = [
      0.0 0.0 0.0 0.0
      0.0 0.0 0.0 0.0
      0.0 0.0 1.0 0.0
      0.0 0.0 0.0 1.0
    ]
  M = kron(Ndn, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Nupdn", ::SiteType"HubHolst", s::Index)
  Nupdn = [
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 1.0
  ]
  M = kron(Nupdn, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Ntot", ::SiteType"HubHolst", s::Index)
  Ntot = [
    0.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 2.0
  ]
  M = kron(Ntot, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Cup", ::SiteType"HubHolst", s::Index)
  Cup = [
    0.0 1.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 1.0
    0.0 0.0 0.0 0.0
  ]
  M = kron(Cup, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Cdagup", ::SiteType"HubHolst", s::Index)
  Cdagup = [
    0.0 0.0 0.0 0.0
    1.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 1.0 0.0
  ]
  M = kron(Cdagup, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Cdn", ::SiteType"HubHolst", s::Index)
  Cdn = [
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 -1.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
  ]
  M = kron(Cdn, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Cupdn", ::SiteType"HubHolst", s::Index)
  Cupdn = [
    0.0  0.0  0.0  -1.0
    0.0  0.0  0.0   0.0
    0.0  0.0  0.0   0.0
    0.0  0.0  0.0   0.0
  ]
  M = kron(Cupdn, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Cdagupdn", ::SiteType"HubHolst", s::Index)
  Cdagupdn = [
    0.0  0.0  0.0   0.0
    0.0  0.0  0.0   0.0
    0.0  0.0  0.0   0.0
   -1.0  0.0  0.0   0.0
  ]
  M = kron(Cdagupdn, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Cdagdn", ::SiteType"HubHolst", s::Index)
  Cdagdn = [
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    1.0 0.0 0.0 0.0
    0.0 -1.0 0.0 0.0
  ]
  M = kron(Cdagdn, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Aup", ::SiteType"HubHolst", s::Index)
  Aup = [
    0.0 1.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 1.0
    0.0 0.0 0.0 0.0
  ]
  M = kron(Aup, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Adagup", ::SiteType"HubHolst", s::Index)
  Adagup = [
    0.0 0.0 0.0 0.0
    1.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 1.0 0.0
  ]
  M = kron(Adagup, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Adn", ::SiteType"HubHolst", s::Index)
  Adn = [
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 -1.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
  ]
  M = kron(Adn, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Adagdn", ::SiteType"HubHolst", s::Index)
  Adagdn = [
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    1.0 0.0 0.0 0.0
    0.0 -1.0 0.0 0.0
  ]
  M = kron(Adagdn, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"F", ::SiteType"HubHolst", s::Index)
  F = [
    1.0 0.0 0.0 0.0
    0.0 -1.0 0.0 0.0
    0.0 0.0 -1.0 0.0
    0.0 0.0 0.0 1.0
  ]
  M = kron(F, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Fup", ::SiteType"HubHolst", s::Index)
  Fup = [
    1.0 0.0 0.0 0.0
    0.0 -1.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 -1.0
  ]
  M = kron(Fup, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Fdn", ::SiteType"HubHolst", s::Index)
  Fdn = [
    1.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 -1.0 0.0
    0.0 0.0 0.0 -1.0
  ]
  M = kron(Fdn, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Sz", ::SiteType"HubHolst", s::Index)
  Sz = [
    0.0 0.0 0.0 0.0
    0.0 0.5 0.0 0.0
    0.0 0.0 -0.5 0.0
    0.0 0.0 0.0 0.0
  ]
  M = kron(Sz, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Sx", ::SiteType"HubHolst", s::Index)
  Sx = [
    0.0 0.0 0.0 0.0
    0.0 0.0 0.5 0.0
    0.0 0.5 0.0 0.0
    0.0 0.0 0.0 0.0
  ]
  M = kron(Sx, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"S+", ::SiteType"HubHolst", s::Index)
  Splus = [
    0.0 0.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
  ]
  M = kron(Splus, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"S-", ::SiteType"HubHolst", s::Index)
  Sminus = [
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 0.0 0.0
  ]
  M = kron(Sminus, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

ITensors.has_fermion_string(::OpName"Cup", ::SiteType"HubHolst") = true
function ITensors.has_fermion_string(on::OpName"c↑", st::SiteType"HubHolst")
  return has_fermion_string(alias(on), st)
end
ITensors.has_fermion_string(::OpName"Cdagup", ::SiteType"HubHolst") = true
function ITensors.has_fermion_string(on::OpName"c†↑", st::SiteType"HubHolst")
  return has_fermion_string(alias(on), st)
end
ITensors.has_fermion_string(::OpName"Cdn", ::SiteType"HubHolst") = true
function ITensors.has_fermion_string(on::OpName"c↓", st::SiteType"HubHolst")
  return has_fermion_string(alias(on), st)
end
ITensors.has_fermion_string(::OpName"Cdagdn", ::SiteType"HubHolst") = true
function ITensors.has_fermion_string(on::OpName"c†↓", st::SiteType"HubHolst")
  return has_fermion_string(alias(on), st)
end

## PHONON OPERATORS ## 

function ITensors.op(::OpName"B", ::SiteType"HubHolst", s::Index)
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons"

  B = Tridiagonal(zeros(floor(Int,dim(s)/4)-1),zeros(floor(Int,dim(s)/4)),sqrt.(collect(1:floor(Int,dim(s)/4)-1)))
  M = kron(Matrix(I, 4, 4), B)
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Bdag", ::SiteType"HubHolst", s::Index)
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons"

  Bdag = sparse(Tridiagonal(sqrt.(collect(1:floor(Int,dim(s)/4)-1)),zeros(floor(Int,dim(s)/4)),zeros(floor(Int,dim(s)/4)-1)))
  M = kron(Matrix(I, 4, 4), Bdag)
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Nb", ::SiteType"HubHolst", s::Index)
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons"

  Bdag = sparse(Tridiagonal(sqrt.(collect(1:floor(Int,dim(s)/4)-1)),zeros(floor(Int,dim(s)/4)),zeros(floor(Int,dim(s)/4)-1)))
  B = Tridiagonal(zeros(floor(Int,dim(s)/4)-1),zeros(floor(Int,dim(s)/4)),sqrt.(collect(1:floor(Int,dim(s)/4)-1)))
  M = kron(Matrix(I, 4, 4), (Bdag * B))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Bdag+B", ::SiteType"HubHolst", s::Index)
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons"

  Bdag = sparse(Tridiagonal(sqrt.(collect(1:floor(Int,dim(s)/4)-1)),zeros(floor(Int,dim(s)/4)),zeros(floor(Int,dim(s)/4)-1)))
  B = Tridiagonal(zeros(floor(Int,dim(s)/4)-1),zeros(floor(Int,dim(s)/4)),sqrt.(collect(1:floor(Int,dim(s)/4)-1)))
  M = kron(Matrix(I, 4, 4), (Bdag + B))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Ntot(Bd+B)", ::SiteType"HubHolst", s::Index)
  # Fermion part 
  Ntot = [
    0.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 2.0
  ] 
  
  # Make B†+B
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons" 

  Bdag = sparse(Tridiagonal(sqrt.(collect(1:floor(Int,dim(s)/4)-1)),zeros(floor(Int,dim(s)/4)),zeros(floor(Int,dim(s)/4)-1))) 
  B = Tridiagonal(zeros(floor(Int,dim(s)/4)-1),zeros(floor(Int,dim(s)/4)),sqrt.(collect(1:floor(Int,dim(s)/4)-1)))
  M = kron(Ntot, (Bdag + B))
  ITensor(M, prime(s), dag(s))
end
using LinearAlgebra
using SparseArrays

PHONON_DOF = 2
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
  conserve_qns=true,
  conserve_sz=conserve_qns,
  conserve_nf=conserve_qns,
  conserve_nfparity=conserve_qns,
  qnname_sz="Sz",
  qnname_nf="Nf",
  qnname_nfparity="NfParity",
)
  if conserve_sz && conserve_nf
    println("Conserving Sz and Nf")
    return [
      QN((qnname_nf, 0, -1), (qnname_sz, 0)) => 1 * PHONON_DOF 
      QN((qnname_nf, 1, -1), (qnname_sz, +1)) => 1 * PHONON_DOF
      QN((qnname_nf, 1, -1), (qnname_sz, -1)) => 1 * PHONON_DOF
      QN((qnname_nf, 2, -1), (qnname_sz, 0)) => 1 * PHONON_DOF
    ]
  elseif conserve_nf
    return [
      QN(qnname_nf, 0, -1) => 1 * PHONON_DOF
      QN(qnname_nf, 1, -1) => 2 * PHONON_DOF
      QN(qnname_nf, 2, -1) => 1 * PHONON_DOF
    ]
  elseif conserve_sz
    return [
      QN((qnname_sz, 0), (qnname_nfparity, 0, -2)) => 1 * PHONON_DOF
      QN((qnname_sz, +1), (qnname_nfparity, 1, -2)) => 1 * PHONON_DOF
      QN((qnname_sz, -1), (qnname_nfparity, 1, -2)) => 1 * PHONON_DOF
      QN((qnname_sz, 0), (qnname_nfparity, 0, -2)) => 1 * PHONON_DOF
    ]
  elseif conserve_nfparity
    return [
      QN(qnname_nfparity, 0, -2) => 1 * PHONON_DOF
      QN(qnname_nfparity, 1, -2) => 2 * PHONON_DOF
      QN(qnname_nfparity, 0, -2) => 1 * PHONON_DOF
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

function ITensors.state(n::StateName{N}, ::SiteType"HubHolst") where {N}
    hubbard_type, phonon_num = split(String(name(n)), ",")
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

    st_ph = zeros(PHONON_DOF)
    st_ph[parse(Int, String(phonon_num)) + 1] = 1

    return kron(st_hubb, st_ph)
end

## ELECTRON OPERATORS ## 

alias(::OpName"c↑") = OpName("Cup")
alias(::OpName"c↓") = OpName("Cdn")
alias(::OpName"c†↑") = OpName("Cdagup")
alias(::OpName"c†↓") = OpName("Cdagdn")
alias(::OpName"n↑") = OpName("Nup")
alias(::OpName"n↓") = OpName("Ndn")
alias(::OpName"n↑↓") = OpName("Nupdn")
alias(::OpName"ntot") = OpName("Ntot")
alias(::OpName"F↑") = OpName("Fup")
alias(::OpName"F↓") = OpName("Fdn")

function ITensors.op(::OpName"Nup", ::SiteType"HubHolst")
  Nup = [
    0.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 1.0
  ]
  return kron(Nup, Matrix(I, PHONON_DOF, PHONON_DOF))
end
function ITensors.op(on::OpName"n↑", st::SiteType"HubHolst")
  return op(alias(on), st)
end

function ITensors.op(::OpName"Ndn", ::SiteType"HubHolst")
  Ndn = [
      0.0 0.0 0.0 0.0
      0.0 0.0 0.0 0.0
      0.0 0.0 1.0 0.0
      0.0 0.0 0.0 1.0
    ]
  return kron(Ndn, Matrix(I, PHONON_DOF, PHONON_DOF))
end
function ITensors.op(on::OpName"n↓", st::SiteType"HubHolst")
  return op(alias(on), st)
end

function ITensors.op(::OpName"Nupdn", ::SiteType"HubHolst")
  Nupdn = [
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 1.0
  ]
  return kron(Nupdn, Matrix(I, PHONON_DOF, PHONON_DOF))
end
function ITensors.op(on::OpName"n↑↓", st::SiteType"HubHolst")
  return op(alias(on), st)
end

function ITensors.op(::OpName"Ntot", ::SiteType"HubHolst")
  Ntot = [
    0.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 2.0
  ]
  return kron(Ntot, Matrix(I, PHONON_DOF, PHONON_DOF))
end
function ITensors.op(on::OpName"ntot", st::SiteType"HubHolst")
  return op(alias(on), st)
end

function ITensors.op(::OpName"Cup", ::SiteType"HubHolst")
  Cup = [
    0.0 1.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 1.0
    0.0 0.0 0.0 0.0
  ]
  return kron(Cup, Matrix(I, PHONON_DOF, PHONON_DOF))
end
function ITensors.op(on::OpName"c↑", st::SiteType"HubHolst")
  return op(alias(on), st)
end

function ITensors.op(::OpName"Cdagup", ::SiteType"HubHolst")
  Cdagup = [
    0.0 0.0 0.0 0.0
    1.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 1.0 0.0
  ]
  return kron(Cdagup, Matrix(I, PHONON_DOF, PHONON_DOF))
end
function ITensors.op(on::OpName"c†↑", st::SiteType"HubHolst")
  return op(alias(on), st)
end

function ITensors.op(::OpName"Cdn", ::SiteType"HubHolst")
  Cdn = [
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 -1.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
  ]
  return kron(Cdn, Matrix(I, PHONON_DOF, PHONON_DOF))
end
function ITensors.op(on::OpName"c↓", st::SiteType"HubHolst")
  return op(alias(on), st)
end

function ITensors.op(::OpName"Cdagdn", ::SiteType"HubHolst")
  Cdagdn = [
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    1.0 0.0 0.0 0.0
    0.0 -1.0 0.0 0.0
  ]
  return kron(Cdagdn, Matrix(I, PHONON_DOF, PHONON_DOF))
end
function ITensors.op(::OpName"c†↓", st::SiteType"HubHolst")
  return op(OpName("Cdagdn"), st)
end

function ITensors.op(::OpName"Aup", ::SiteType"HubHolst")
  Aup = [
    0.0 1.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 1.0
    0.0 0.0 0.0 0.0
  ]
  return kron(Aup, Matrix(I, PHONON_DOF, PHONON_DOF))
end
function ITensors.op(on::OpName"a↑", st::SiteType"HubHolst")
  return op(alias(on), st)
end

function ITensors.op(::OpName"Adagup", ::SiteType"HubHolst")
  Adagup = [
    0.0 0.0 0.0 0.0
    1.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 1.0 0.0
  ]
  return kron(Adagup, Matrix(I, PHONON_DOF, PHONON_DOF))
end
function ITensors.op(on::OpName"a†↑", st::SiteType"HubHolst")
  return op(alias(on), st)
end

function ITensors.op(::OpName"Adn", ::SiteType"HubHolst")
  Adn = [
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 -1.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
  ]
  return kron(Adn, Matrix(I, PHONON_DOF, PHONON_DOF))
end
function ITensors.op(on::OpName"a↓", st::SiteType"HubHolst")
  return op(alias(on), st)
end

function ITensors.op(::OpName"Adagdn", ::SiteType"HubHolst")
  Adagdn = [
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    1.0 0.0 0.0 0.0
    0.0 -1.0 0.0 0.0
  ]
  return kron(Adagdn, Matrix(I, PHONON_DOF, PHONON_DOF))
end
function ITensors.op(::OpName"a†↓", st::SiteType"HubHolst")
  return op(OpName("Cdagdn"), st)
end


function ITensors.op(::OpName"F", ::SiteType"HubHolst")
  F = [
    1.0 0.0 0.0 0.0
    0.0 -1.0 0.0 0.0
    0.0 0.0 -1.0 0.0
    0.0 0.0 0.0 1.0
  ]
  return kron(F, Matrix(I, PHONON_DOF, PHONON_DOF))
end

function ITensors.op(::OpName"Fup", ::SiteType"HubHolst")
  Fup = return [
    1.0 0.0 0.0 0.0
    0.0 -1.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 -1.0
  ]
  return kron(Fup, Matrix(I, PHONON_DOF, PHONON_DOF))
end
function ITensors.op(::OpName"F↑", st::SiteType"HubHolst")
  return op(OpName("Fup"), st)
end

function ITensors.op(::OpName"Fdn", ::SiteType"HubHolst")
  Fdn = [
    1.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 -1.0 0.0
    0.0 0.0 0.0 -1.0
  ]
  return kron(Fdn, Matrix(I, PHONON_DOF, PHONON_DOF))
end
function ITensors.op(::OpName"F↓", st::SiteType"HubHolst")
  return op(OpName("Fdn"), st)
end

function ITensors.op(::OpName"Sz", ::SiteType"HubHolst")
  #Op[s' => 2, s => 2] = +0.5
  #return Op[s' => 3, s => 3] = -0.5
  Sz = [
    0.0 0.0 0.0 0.0
    0.0 0.5 0.0 0.0
    0.0 0.0 -0.5 0.0
    0.0 0.0 0.0 0.0
  ]
  return kron(Sz, Matrix(I, PHONON_DOF, PHONON_DOF))
end

function ITensors.op(::OpName"Sᶻ", st::SiteType"HubHolst")
  return op(OpName("Sz"), st)
end

function ITensors.op(::OpName"Sx", ::SiteType"HubHolst")
  Sx = [
    0.0 0.0 0.0 0.0
    0.0 0.0 0.5 0.0
    0.0 0.5 0.0 0.0
    0.0 0.0 0.0 0.0
  ]
  return kron(Sx, Matrix(I, PHONON_DOF, PHONON_DOF))
end

function ITensors.op(::OpName"Sˣ", st::SiteType"HubHolst")
  return op(OpName("Sx"), st)
end

function ITensors.op(::OpName"S+", ::SiteType"HubHolst")
  Splus = [
    0.0 0.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
  ]
  return kron(Splus, Matrix(I, PHONON_DOF, PHONON_DOF))
end

function ITensors.op(::OpName"S⁺", st::SiteType"HubHolst")
  return op(OpName("S+"), st)
end
function ITensors.op(::OpName"Sp", st::SiteType"HubHolst")
  return op(OpName("S+"), st)
end
function ITensors.op(::OpName"Splus", st::SiteType"HubHolst")
  return op(OpName("S+"), st)
end

function ITensors.op(::OpName"S-", ::SiteType"HubHolst")
  Sminus = [
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 0.0 0.0
  ]
  return kron(Sminus, Matrix(I, PHONON_DOF, PHONON_DOF))
end

function ITensors.op(::OpName"S⁻", st::SiteType"HubHolst")
  return op(OpName("S-"), st)
end
function ITensors.op(::OpName"Sm", st::SiteType"HubHolst")
  return op(OpName("S-"), st)
end
function ITensors.op(::OpName"Sminus", st::SiteType"HubHolst")
  return op(OpName("S-"), st)
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

function ITensors.op(::OpName"I", ::SiteType"HubHolst")
  return Matrix(I, 4*PHONON_DOF, 4*PHONON_DOF) # Total state space is 4 (from Hubbard) * PHONON_DOF (# phonon dof)
end

function ITensors.op(::OpName"B", ::SiteType"HubHolst")
  if PHONON_DOF > 1
    B = Tridiagonal(zeros(PHONON_DOF-1),zeros(PHONON_DOF),sqrt.(collect(1:PHONON_DOF-1)))
    return kron(Matrix(I, 4, 4), B)
  end
  return Matrix(I, 4*PHONON_DOF, 4*PHONON_DOF)
end

function ITensors.op(::OpName"Bdag", ::SiteType"HubHolst")
  if PHONON_DOF > 1
    Bdag = sparse(Tridiagonal(sqrt.(collect(1:PHONON_DOF-1)),zeros(PHONON_DOF),zeros(PHONON_DOF-1)))
    return kron(Matrix(I, 4, 4), Bdag)
  end
  return Matrix(I, 4*PHONON_DOF, 4*PHONON_DOF)
end

function ITensors.op(::OpName"Nb", ::SiteType"HubHolst")
  if PHONON_DOF > 1
    Bdag = sparse(Tridiagonal(sqrt.(collect(1:PHONON_DOF-1)),zeros(PHONON_DOF),zeros(PHONON_DOF-1)))
    B = Tridiagonal(zeros(PHONON_DOF-1),zeros(PHONON_DOF),sqrt.(collect(1:PHONON_DOF-1)))
    return kron(Matrix(I, 4, 4), (Bdag * B))
  end
  return 0#Matrix(I, 4*PHONON_DOF, 4*PHONON_DOF)
end

function ITensors.op(::OpName"Bdag+B", st::SiteType"HubHolst")
  Id = Matrix(I, 4, 4) # Hubbard model has d=4 for fermion sites
  if PHONON_DOF > 1
    Bdag = sparse(Tridiagonal(sqrt.(collect(1:PHONON_DOF-1)),zeros(PHONON_DOF),zeros(PHONON_DOF-1)))
    B = Tridiagonal(zeros(PHONON_DOF-1),zeros(PHONON_DOF),sqrt.(collect(1:PHONON_DOF-1)))
    return kron(Matrix(I, 4, 4), (Bdag + B))
  end
  return Matrix(I, 4*PHONON_DOF, 4*PHONON_DOF)
end

function ITensors.op(::OpName"Ntot(Bd+B)", st::SiteType"HubHolst")
  Nf = [
    0.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 2.0
  ] # Hubbard model has d=4 for fermion sites
  
  # Make B† 
  if PHONON_DOF > 1
    Bdag = sparse(Tridiagonal(sqrt.(collect(1:PHONON_DOF-1)),zeros(PHONON_DOF),zeros(PHONON_DOF-1))) 
    B = Tridiagonal(zeros(PHONON_DOF-1),zeros(PHONON_DOF),sqrt.(collect(1:PHONON_DOF-1)))
    return kron(Matrix(I, 4, 4), (Bdag + B))
  end
  return 0#Matrix(I, 4*PHONON_DOF, 4*PHONON_DOF)
end
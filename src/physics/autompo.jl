#
# Optimizations:
#  - replace leftmap, rightmap with sorted vectors
# 

###########################
# SiteOp                  # 
###########################

struct SiteOp{N}
  name::String
  site::NTuple{N,Int}
  params::NamedTuple
end
# Change NamedTuple() to (;) when we drop older Julia versions
SiteOp(name::String, site::Tuple) = SiteOp(name, site, NamedTuple())
SiteOp(name::String, site::Int...) = SiteOp(name, site)
function SiteOp(name::String, site_params::Union{Int,NamedTuple}...)
  return SiteOp(name, Base.front(site_params), last(site_params))
end
SiteOp(name::String, params::NamedTuple, site::Tuple) = SiteOp(name, site, params)
SiteOp(name::String, params::NamedTuple, site::Int...) = SiteOp(name, site, params)

convert(::Type{SiteOp}, op::Pair{String,Int}) = SiteOp(first(op), last(op))

name(s::SiteOp) = s.name
site(s::SiteOp) = only(s.site)
sites(s::SiteOp) = s.site
params(s::SiteOp) = s.params

site_or_sites(s::SiteOp{1}) = site(s)
site_or_sites(s::SiteOp) = sites(s)

string_site_or_sites(s::SiteOp{1}) = string(site(s))
string_site_or_sites(s::SiteOp) = string(sites(s))[2:(end - 1)]

show(io::IO, s::SiteOp) = print(io, "\"$(name(s))\"($(string_site_or_sites(s)))")

(s1::SiteOp == s2::SiteOp) = (s1.site == s2.site && s1.name == s2.name)

function isless(s1::SiteOp, s2::SiteOp)
  if site(s1) != site(s2)
    return site(s1) < site(s2)
  end
  return name(s1) < name(s2)
end

###########################
# OpTerm                  # 
###########################

const OpTerm = Vector{SiteOp}

function (o1::OpTerm == o2::OpTerm)
  (length(o1) == length(o2)) || return false
  @inbounds for n in 1:length(o1)
    (o1[n] != o2[n]) && return false
  end
  return true
end

function isless(o1::OpTerm, o2::OpTerm)
  if length(o1) != length(o2)
    return length(o1) < length(o2)
  end
  for n in 1:length(o1)
    if o1[n] != o2[n]
      return (o1[n] < o2[n])
    end
  end
  return false
end

mult(t1::OpTerm, t2::OpTerm) = isempty(t2) ? t1 : vcat(t1, t2)

function isfermionic(t::OpTerm, sites)::Bool
  p = +1
  for op in t
    if has_fermion_string(name(op), sites[site(op)])
      p *= -1
    end
  end
  return (p == -1)
end

###########################
# MPOTerm                 # 
###########################

mutable struct MPOTerm
  coef::ComplexF64
  ops::OpTerm
end
coef(op::MPOTerm) = op.coef
ops(op::MPOTerm) = op.ops

copy(t::MPOTerm) = MPOTerm(coef(t), copy(ops(t)))

function (t1::MPOTerm == t2::MPOTerm)
  return coef(t1) ≈ coef(t2) && ops(t1) == ops(t2)
end

function isless(t1::MPOTerm, t2::MPOTerm)
  if ops(t1) == ops(t2)
    if coef(t1) ≈ coef(t2)
      return false
    else
      ct1 = coef(t1)
      ct2 = coef(t2)
      #"lexicographic" ordering on  complex numbers
      return real(ct1) < real(ct2) || (real(ct1) ≈ real(ct2) && imag(ct1) < imag(ct2))
    end
  end
  return ops(t1) < ops(t2)
end

function MPOTerm(c::Number, op1::String, ops_rest...)
  ops = (op1, ops_rest...)
  starts = findall(x -> x isa String, ops)
  N = length(starts)
  vop = OpTerm(undef, N)
  for n in 1:N
    start = starts[n]
    stop = (n == N) ? lastindex(ops) : (starts[n + 1] - 1)
    vop[n] = SiteOp(ops[start:stop]...)
  end
  return MPOTerm(c, vop)
end

function MPOTerm(op1::String, ops...)
  return MPOTerm(one(Float64), op1, ops...)
end

function MPOTerm(ops::Vector{Pair{String,Int}})
  return MPOTerm(Iterators.flatten(ops)...)
end

function Base.show(io::IO, op::MPOTerm)
  c = coef(op)
  if iszero(imag(c))
    print(io, "$(real(c)) ")
  elseif iszero(real(c))
    print(io, "$(imag(c))im ")
  else
    print(io, "($c) ")
  end
  for o in ops(op)
    print(io, "\"$(name(o))\"($(string_site_or_sites(o))) ")
  end
end

############################
## OpSum                 #
############################

"""
An `OpSum` represents a sum of operator
terms.

Often it is used to create matrix
product operator (`MPO`) approximation
of the sum of the terms in the `OpSum` oject.
Each term is a product of local operators
specified by names such as `"Sz"` or `"N"`,
times an optional coefficient which
can be real or complex.

Which local operator names are available
is determined by the function `op`
associated with the `TagType` defined by
special Index tags, such as `"S=1/2"`, `"S=1"`,
`"Fermion"`, and `"Electron"`.
"""
mutable struct OpSum
  data::Vector{MPOTerm}
  OpSum(terms::Vector{MPOTerm}) = new(terms)
end

length(os::OpSum) = length(data(os))
getindex(os::OpSum, I...) = data(os)[I...]

const AutoMPO = OpSum

"""
    OpSum()
    
Construct an empty `OpSum`.
"""
OpSum() = OpSum(Vector{MPOTerm}())

data(ampo::OpSum) = ampo.data
setdata!(ampo::OpSum, ndata) = (ampo.data = ndata)

push!(ampo::OpSum, term) = push!(data(ampo), term)

Base.:(==)(ampo1::OpSum, ampo2::OpSum) = data(ampo1) == data(ampo2)

Base.copy(ampo::OpSum) = OpSum(copy(data(ampo)))

function Base.deepcopy(ampo::OpSum)
  return OpSum(map(copy, data(ampo)))
end

Base.size(ampo::OpSum) = size(data(ampo))

"""
    add!(ampo::OpSum,
         op1::String, i1::Int)

    add!(ampo::OpSum,
         coef::Number,
         op1::String, i1::Int)

    add!(ampo::OpSum,
         op1::String, i1::Int,
         op2::String, i2::Int,
         ops...)

    add!(ampo::OpSum,
         coef::Number,
         op1::String, i1::Int,
         op2::String, i2::Int,
         ops...)

    +(ampo:OpSum, term::Tuple)

Add a single- or multi-site operator 
term to the OpSum `ampo`. Each operator
is specified by a name (String) and a
site number (Int). The second version
accepts a real or complex coefficient.

The `+` operator version of this function
accepts a tuple with entries either
(String,Int,String,Int,...) or
(Number,String,Int,String,Int,...)
where these tuple values are the same
as valid inputs to the `add!` function.
For inputting a very large number of
terms (tuples) to an OpSum, consider
using the broadcasted operator `.+=`
which avoids reallocating the OpSum
after each addition.

# Examples
```julia
ampo = OpSum()

add!(ampo,"Sz",2,"Sz",3)

ampo += ("Sz",3,"Sz",4)

ampo += (0.5,"S+",4,"S-",5)

ampo .+= (0.5,"S+",5,"S-",6)
```
"""
add!(os::OpSum, t::MPOTerm) = push!(os, t)

add!(os::OpSum, args...) = add!(os, MPOTerm(args...))

"""
    subtract!(ampo::OpSum,
              op1::String, i1::Int,
              op2::String, i2::Int,
              ops...)

    subtract!(ampo::OpSum,
              coef::Number,
              op1::String, i1::Int,
              op2::String, i2::Int,
              ops...)

Subtract a multi-site operator term
from the OpSum `ampo`. Each operator
is specified by a name (String) and a
site number (Int). The second version
accepts a real or complex coefficient.
"""
subtract!(os::OpSum, args...) = add!(os, -MPOTerm(args...))

-(t::MPOTerm) = MPOTerm(-coef(t), ops(t))

function (ampo::OpSum + term::Tuple)
  ampo_plus_term = copy(ampo)
  add!(ampo_plus_term, term...)
  return ampo_plus_term
end

function (ampo::OpSum + term::Vector{Pair{String,Int64}})
  ampo_plus_term = copy(ampo)
  add!(ampo_plus_term, term)
  return ampo_plus_term
end

function (ampo::OpSum - term::Tuple)
  ampo_plus_term = copy(ampo)
  subtract!(ampo_plus_term, term...)
  return ampo_plus_term
end

#
# ampo .+= ("Sz",1) syntax using broadcasting
#

struct OpSumStyle <: Broadcast.BroadcastStyle end
Base.BroadcastStyle(::Type{<:OpSum}) = OpSumStyle()

struct OpSumAddTermStyle <: Broadcast.BroadcastStyle end

Base.broadcastable(ampo::OpSum) = ampo

Base.BroadcastStyle(::OpSumStyle, ::Broadcast.Style{Tuple}) = OpSumAddTermStyle()

Broadcast.instantiate(bc::Broadcast.Broadcasted{OpSumAddTermStyle}) = bc

function Base.copyto!(ampo, bc::Broadcast.Broadcasted{OpSumAddTermStyle,<:Any,typeof(+)})
  add!(ampo, bc.args[2]...)
  return ampo
end

#
# ampo .-= ("Sz",1) syntax using broadcasting
#

function Base.copyto!(ampo, bc::Broadcast.Broadcasted{OpSumAddTermStyle,<:Any,typeof(-)})
  subtract!(ampo, bc.args[2]...)
  return ampo
end

function Base.show(io::IO, ampo::OpSum)
  println(io, "OpSum:")
  for term in data(ampo)
    println(io, "  $term")
  end
end

##################################
# MatElem (simple sparse matrix) #
##################################

struct MatElem{T}
  row::Int
  col::Int
  val::T
end

#function Base.show(io::IO,m::MatElem)
#  print(io,"($(m.row),$(m.col),$(m.val))")
#end

function toMatrix(els::Vector{MatElem{T}})::Matrix{T} where {T}
  nr = 0
  nc = 0
  for el in els
    nr = max(nr, el.row)
    nc = max(nc, el.col)
  end
  M = zeros(T, nr, nc)
  for el in els
    M[el.row, el.col] = el.val
  end
  return M
end

function Base.:(==)(m1::MatElem{T}, m2::MatElem{T})::Bool where {T}
  return (m1.row == m2.row && m1.col == m2.col && m1.val == m2.val)
end

function Base.isless(m1::MatElem{T}, m2::MatElem{T})::Bool where {T}
  if m1.row != m2.row
    return m1.row < m2.row
  elseif m1.col != m2.col
    return m1.col < m2.col
  end
  return m1.val < m2.val
end

struct QNMatElem{T}
  rowqn::QN
  colqn::QN
  row::Int
  col::Int
  val::T
end

function Base.:(==)(m1::QNMatElem{T}, m2::QNMatElem{T})::Bool where {T}
  return (
    m1.row == m2.row &&
    m1.col == m2.col &&
    m1.val == m2.val &&
    m1.rowqn == m2.rowqn &&
    m1.colqn == m2.colqn
  )
end

function Base.isless(m1::QNMatElem{T}, m2::QNMatElem{T})::Bool where {T}
  if m1.rowqn != m2.rowqn
    return m1.rowqn < m2.rowqn
  elseif m1.colqn != m2.colqn
    return m1.colqn < m2.colqn
  elseif m1.row != m2.row
    return m1.row < m2.row
  elseif m1.col != m2.col
    return m1.col < m2.col
  end
  return m1.val < m2.val
end

isempty(op_qn::Pair{OpTerm,QN}) = isempty(op_qn.first)

# the key type is OpTerm for the dense case
# and is Pair{OpTerm,QN} for the QN conserving case
function posInLink!(linkmap::Dict{K,Int}, k::K)::Int where {K}
  isempty(k) && return -1
  pos = get(linkmap, k, -1)
  if pos == -1
    pos = length(linkmap) + 1
    linkmap[k] = pos
  end
  return pos
end

function determineValType(terms::Vector{MPOTerm})
  for t in terms
    (!isreal(coef(t))) && return ComplexF64
  end
  return Float64
end

function computeSiteProd(sites, ops::OpTerm)::ITensor
  i = site(ops[1])
  T = op(sites[i], ops[1].name)
  for j in 2:length(ops)
    (site(ops[j]) != i) && error("Mismatch of site number in computeSiteProd")
    opj = op(sites[i], ops[j].name)
    T = product(T, opj)
  end
  return T
end

function remove_dups!(v::Vector{T}) where {T}
  N = length(v)
  (N == 0) && return nothing
  sort!(v)
  n = 1
  u = 2
  while u <= N
    while u < N && v[u] == v[n]
      u += 1
    end
    if v[u] != v[n]
      v[n + 1] = v[u]
      n += 1
    end
    u += 1
  end
  resize!(v, n)
  return nothing
end #remove_dups!

function svdMPO(ampo::OpSum, sites; kwargs...)::MPO
  mindim::Int = get(kwargs, :mindim, 1)
  maxdim::Int = get(kwargs, :maxdim, 10000)
  cutoff::Float64 = get(kwargs, :cutoff, 1E-15)

  N = length(sites)

  ValType = determineValType(data(ampo))

  Vs = [Matrix{ValType}(undef, 1, 1) for n in 1:N]
  tempMPO = [MatElem{MPOTerm}[] for n in 1:N]

  crosses_bond(t::MPOTerm, n::Int) = (site(ops(t)[1]) <= n <= site(ops(t)[end]))

  rightmap = Dict{OpTerm,Int}()
  next_rightmap = Dict{OpTerm,Int}()

  for n in 1:N
    leftbond_coefs = MatElem{ValType}[]

    leftmap = Dict{OpTerm,Int}()
    for term in data(ampo)
      crosses_bond(term, n) || continue

      left::OpTerm = filter(t -> (site(t) < n), ops(term))
      onsite::OpTerm = filter(t -> (site(t) == n), ops(term))
      right::OpTerm = filter(t -> (site(t) > n), ops(term))

      bond_row = -1
      bond_col = -1
      if !isempty(left)
        bond_row = posInLink!(leftmap, left)
        bond_col = posInLink!(rightmap, mult(onsite, right))
        bond_coef = convert(ValType, coef(term))
        push!(leftbond_coefs, MatElem(bond_row, bond_col, bond_coef))
      end

      A_row = bond_col
      A_col = posInLink!(next_rightmap, right)
      site_coef = 1.0 + 0.0im
      if A_row == -1
        site_coef = coef(term)
      end
      if isempty(onsite)
        if !using_auto_fermion() && isfermionic(right, sites)
          push!(onsite, SiteOp("F", n))
        else
          push!(onsite, SiteOp("Id", n))
        end
      end
      el = MatElem(A_row, A_col, MPOTerm(site_coef, onsite))
      push!(tempMPO[n], el)
    end
    rightmap = next_rightmap
    next_rightmap = Dict{OpTerm,Int}()

    remove_dups!(tempMPO[n])

    if n > 1 && !isempty(leftbond_coefs)
      M = toMatrix(leftbond_coefs)
      U, S, V = svd(M)
      P = S .^ 2
      truncate!(P; maxdim=maxdim, cutoff=cutoff, mindim=mindim)
      tdim = length(P)
      nc = size(M, 2)
      Vs[n - 1] = Matrix{ValType}(V[1:nc, 1:tdim])
    end
  end

  llinks = Vector{Index{Int}}(undef, N + 1)
  llinks[1] = Index(2, "Link,l=0")

  H = MPO(sites)

  for n in 1:N
    VL = Matrix{ValType}(undef, 1, 1)
    if n > 1
      VL = Vs[n - 1]
    end
    VR = Vs[n]
    tdim = size(VR, 2)

    llinks[n + 1] = Index(2 + tdim, "Link,l=$n")

    ll = llinks[n]
    rl = llinks[n + 1]

    H[n] = ITensor()

    for el in tempMPO[n]
      A_row = el.row
      A_col = el.col
      t = el.val
      (abs(coef(t)) > eps()) || continue

      M = zeros(ValType, dim(ll), dim(rl))

      ct = convert(ValType, coef(t))
      if A_row == -1 && A_col == -1 #onsite term
        M[end, 1] += ct
      elseif A_row == -1 #term starting on site n
        for c in 1:size(VR, 2)
          z = ct * VR[A_col, c]
          M[end, 1 + c] += z
        end
      elseif A_col == -1 #term ending on site n
        for r in 1:size(VL, 2)
          z = ct * conj(VL[A_row, r])
          M[1 + r, 1] += z
        end
      else
        for r in 1:size(VL, 2), c in 1:size(VR, 2)
          z = ct * conj(VL[A_row, r]) * VR[A_col, c]
          M[1 + r, 1 + c] += z
        end
      end

      T = itensor(M, ll, rl)
      H[n] += T * computeSiteProd(sites, ops(t))
    end

    #
    # Special handling of starting and 
    # ending identity operators:
    #
    idM = zeros(ValType, dim(ll), dim(rl))
    idM[1, 1] = 1.0
    idM[end, end] = 1.0
    T = itensor(idM, ll, rl)
    H[n] += T * computeSiteProd(sites, SiteOp[SiteOp("Id", n)])
  end

  L = ITensor(llinks[1])
  L[end] = 1.0

  R = ITensor(llinks[N + 1])
  R[1] = 1.0

  H[1] *= L
  H[N] *= R

  return H
end #svdMPO

function qn_svdMPO(ampo::OpSum, sites; kwargs...)::MPO
  mindim::Int = get(kwargs, :mindim, 1)
  maxdim::Int = get(kwargs, :maxdim, 10000)
  cutoff::Float64 = get(kwargs, :cutoff, 1E-15)

  N = length(sites)

  ValType = determineValType(data(ampo))

  Vs = [Dict{QN,Matrix{ValType}}() for n in 1:(N + 1)]
  tempMPO = [QNMatElem{MPOTerm}[] for n in 1:N]

  crosses_bond(t::MPOTerm, n::Int) = (site(ops(t)[1]) <= n <= site(ops(t)[end]))

  rightmap = Dict{Pair{OpTerm,QN},Int}()
  next_rightmap = Dict{Pair{OpTerm,QN},Int}()

  # A cache of the ITensor operators on a certain site
  # of a certain type
  op_cache = Dict{Pair{String,Int},ITensor}()

  for n in 1:N
    leftbond_coefs = Dict{QN,Vector{MatElem{ValType}}}()

    leftmap = Dict{Pair{OpTerm,QN},Int}()
    for term in data(ampo)
      crosses_bond(term, n) || continue

      left::OpTerm = filter(t -> (site(t) < n), ops(term))
      onsite::OpTerm = filter(t -> (site(t) == n), ops(term))
      right::OpTerm = filter(t -> (site(t) > n), ops(term))

      function calcQN(term::OpTerm)
        q = QN()
        for st in term
          op_tensor = get(op_cache, name(st) => site(st), nothing)
          if op_tensor === nothing
            op_tensor = op(sites[site(st)], name(st))
            op_cache[name(st) => site(st)] = op_tensor
          end
          q -= flux(op_tensor)
        end
        return q
      end
      lqn = calcQN(left)
      sqn = calcQN(onsite)

      bond_row = -1
      bond_col = -1
      if !isempty(left)
        bond_row = posInLink!(leftmap, left => lqn)
        bond_col = posInLink!(rightmap, mult(onsite, right) => lqn)
        bond_coef = convert(ValType, coef(term))
        q_leftbond_coefs = get!(leftbond_coefs, lqn, MatElem{ValType}[])
        push!(q_leftbond_coefs, MatElem(bond_row, bond_col, bond_coef))
      end

      rqn = sqn + lqn
      A_row = bond_col
      A_col = posInLink!(next_rightmap, right => rqn)
      site_coef = 1.0 + 0.0im
      if A_row == -1
        site_coef = coef(term)
      end
      if isempty(onsite)
        if !using_auto_fermion() && isfermionic(right, sites)
          push!(onsite, SiteOp("F", n))
        else
          push!(onsite, SiteOp("Id", n))
        end
      end
      el = QNMatElem(lqn, rqn, A_row, A_col, MPOTerm(site_coef, onsite))
      push!(tempMPO[n], el)
    end
    rightmap = next_rightmap
    next_rightmap = Dict{Pair{OpTerm,QN},Int}()

    remove_dups!(tempMPO[n])

    if n > 1 && !isempty(leftbond_coefs)
      for (q, mat) in leftbond_coefs
        M = toMatrix(mat)
        U, S, V = svd(M)
        P = S .^ 2
        truncate!(P; maxdim=maxdim, cutoff=cutoff, mindim=mindim)
        tdim = length(P)
        nc = size(M, 2)
        Vs[n][q] = Matrix{ValType}(V[1:nc, 1:tdim])
      end
    end
  end

  #
  # Make MPO link indices
  #
  d0 = 2
  llinks = Vector{QNIndex}(undef, N + 1)
  # Set dir=In for fermionic ordering, avoid arrow sign
  # <fermions>:
  linkdir = using_auto_fermion() ? In : Out
  llinks[1] = Index(QN() => d0; tags="Link,l=0", dir=linkdir)
  for n in 1:N
    qi = Vector{Pair{QN,Int}}()
    if !haskey(Vs[n + 1], QN())
      # Make sure QN=zero is first in list of sectors
      push!(qi, QN() => d0)
    end
    for (q, Vq) in Vs[n + 1]
      cols = size(Vq, 2)
      if q == QN()
        # Make sure QN=zero is first in list of sectors
        insert!(qi, 1, q => d0 + cols)
      else
        if using_auto_fermion() # <fermions>
          push!(qi, (-q) => cols)
        else
          push!(qi, q => cols)
        end
      end
    end
    # Set dir=In for fermionic ordering, avoid arrow sign
    # <fermions>:
    llinks[n + 1] = Index(qi...; tags="Link,l=$n", dir=linkdir)
  end

  H = MPO(N)

  # Constants which define MPO start/end scheme
  startState = 2
  endState = 1

  for n in 1:N
    finalMPO = Dict{Tuple{QN,OpTerm},Matrix{ValType}}()

    ll = llinks[n]
    rl = llinks[n + 1]

    function defaultMat(ll, rl, lqn, rqn)
      #ldim = qnblockdim(ll,lqn)
      #rdim = qnblockdim(rl,rqn)
      ldim = blockdim(ll, lqn)
      rdim = blockdim(rl, rqn)
      return zeros(ValType, ldim, rdim)
    end

    idTerm = [SiteOp("Id", n)]
    finalMPO[(QN(), idTerm)] = defaultMat(ll, rl, QN(), QN())
    idM = finalMPO[(QN(), idTerm)]
    idM[1, 1] = 1.0
    idM[2, 2] = 1.0

    for el in tempMPO[n]
      t = el.val
      (abs(coef(t)) > eps()) || continue
      A_row = el.row
      A_col = el.col

      M = get!(finalMPO, (el.rowqn, ops(t)), defaultMat(ll, rl, el.rowqn, el.colqn))

      # rowShift and colShift account for
      # special entries in the zero-QN sector
      # of the MPO
      rowShift = (el.rowqn == QN()) ? 2 : 0
      colShift = (el.colqn == QN()) ? 2 : 0

      ct = convert(ValType, coef(t))
      if A_row == -1 && A_col == -1 #onsite term
        M[startState, endState] += ct
      elseif A_row == -1 #term starting on site n
        VR = Vs[n + 1][el.colqn]
        for c in 1:size(VR, 2)
          z = ct * VR[A_col, c]
          M[startState, colShift + c] += z
        end
      elseif A_col == -1 #term ending on site n
        VL = Vs[n][el.rowqn]
        for r in 1:size(VL, 2)
          z = ct * conj(VL[A_row, r])
          M[rowShift + r, endState] += z
        end
      else
        VL = Vs[n][el.rowqn]
        VR = Vs[n + 1][el.colqn]
        for r in 1:size(VL, 2), c in 1:size(VR, 2)
          z = ct * conj(VL[A_row, r]) * VR[A_col, c]
          M[rowShift + r, colShift + c] += z
        end
      end
    end

    s = sites[n]
    H[n] = ITensor()
    for (q_op, M) in finalMPO
      op_prod = q_op[2]
      Op = computeSiteProd(sites, op_prod)

      rq = q_op[1]
      sq = flux(Op)
      cq = rq - sq

      if using_auto_fermion()
        # <fermions>:
        # MPO is defined with Index order
        # of (rl,s[n]',s[n],cl) where rl = row link, cl = col link
        # so compute sign that would result by permuting cl from
        # second position to last position:
        if fparity(sq) == 1 && fparity(cq) == 1
          Op .*= -1
        end
      end

      rn = qnblocknum(ll, rq)
      cn = qnblocknum(rl, cq)

      #TODO: wrap following 3 lines into a function
      _block = Block(rn, cn)
      T = BlockSparseTensor(ValType, [_block], (dag(ll), rl))
      #blockview(T, _block) .= M
      T[_block] .= M

      IT = itensor(T)
      H[n] += IT * Op
    end
  end

  L = ITensor(llinks[1])
  L[startState] = 1.0

  R = ITensor(dag(llinks[N + 1]))
  R[endState] = 1.0

  H[1] *= L
  H[N] *= R

  return H
end #qn_svdMPO

function sorteachterm!(ampo::OpSum, sites)
  ampo = copy(ampo)
  isless_site(o1::SiteOp, o2::SiteOp) = site(o1) < site(o2)
  N = length(sites)
  for t in data(ampo)
    Nt = length(t.ops)
    prevsite = N + 1 #keep track of whether we are switching
    #to a new site to make sure F string
    #is only placed at most once for each site

    # Sort operators in t by site order,
    # and keep the permutation used, perm, for analysis below
    perm = Vector{Int}(undef, Nt)
    sortperm!(perm, t.ops; alg=InsertionSort, lt=isless_site)

    t.ops = t.ops[perm]

    # Identify fermionic operators,
    # zeroing perm for bosonic operators,
    # and inserting string "F" operators
    parity = +1
    for n in Nt:-1:1
      currsite = site(t.ops[n])
      fermionic = has_fermion_string(name(t.ops[n]), sites[site(t.ops[n])])
      if !using_auto_fermion() && (parity == -1) && (currsite < prevsite)
        # Put local piece of Jordan-Wigner string emanating
        # from fermionic operators to the right
        # (Remaining F operators will be put in by svdMPO)
        t.ops[n] = SiteOp("$(name(t.ops[n]))*F", site(t.ops[n]))
      end
      prevsite = currsite

      if fermionic
        parity = -parity
      else
        # Ignore bosonic operators in perm
        # by zeroing corresponding entries
        perm[n] = 0
      end
    end
    if parity == -1
      error("Parity-odd fermionic terms not yet supported by AutoMPO")
    end

    # Keep only fermionic op positions (non-zero entries)
    filter!(!iszero, perm)
    # and account for anti-commuting, fermionic operators 
    # during above sort; put resulting sign into coef
    t.coef *= parity_sign(perm)
  end
  return ampo
end

function sortmergeterms!(ampo::OpSum)
  sort!(data(ampo))

  # Merge (add) terms with same operators
  da = data(ampo)
  ndata = MPOTerm[]
  last_term = copy(da[1])
  for n in 2:length(da)
    if ops(da[n]) == ops(last_term)
      last_term.coef += coef(da[n])
    else
      push!(ndata, last_term)
      last_term = copy(da[n])
    end
  end
  push!(ndata, last_term)

  setdata!(ampo, ndata)
  return ampo
end

"""
    MPO(ampo::OpSum,sites::Vector{<:Index};kwargs...)
       
Convert an OpSum object `ampo` to an
MPO, with indices given by `sites`. The
resulting MPO will have the indices
`sites[1], sites[1]', sites[2], sites[2]'`
etc. The conversion is done by an algorithm
that compresses the MPO resulting from adding
the OpSum terms together, often achieving
the minimum possible bond dimension.

# Examples
```julia
ampo = OpSum()
ampo += ("Sz",1,"Sz",2)
ampo += ("Sz",2,"Sz",3)
ampo += ("Sz",3,"Sz",4)

sites = siteinds("S=1/2",4)
H = MPO(ampo,sites)
```
"""
function MPO(ampo::OpSum, sites::Vector{<:Index}; kwargs...)::MPO
  length(data(ampo)) == 0 && error("OpSum has no terms")

  ampo = deepcopy(ampo)
  sorteachterm!(ampo, sites)
  sortmergeterms!(ampo)

  if hasqns(sites[1])
    return qn_svdMPO(ampo, sites; kwargs...)
  end
  return svdMPO(ampo, sites; kwargs...)
end

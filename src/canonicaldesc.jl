struct CanonicalDescription{T<:AbstractFloat,M<:AbstractMatrix}
  Q::M
  q₁::Vector{T}
  q₂::T
  A::M
  b::Vector{T}
  C̃::M
  c̃::Vector{T}
  lb::Vector{T}
  ub::Vector{T}
  vars::Vector{Symbol}
  name::String

  function (::Type{CanonicalDescription})(qp::MPSDescription{T,M}) where {T<:AbstractFloat,M<:AbstractMatrix}
    c1ind   = isfinite.(qp.c₁)
    c2ind   = isfinite.(qp.c₂)
    m₁, m₂  = sum(c1ind), sum(c2ind)
    m       = m₁ + m₂
    n       = length(qp.vars)
    C̃       = zeros(T, m, n)
    c̃       = zeros(T, m)
    # First, left hand side
    C̃[1:m₁,:]  = -qp.C[c1ind,:]
    c̃[1:m₁]    = -qp.c₁[c1ind]
    # Then, right hand side
    C̃[(m₁+1):end,:]  = qp.C[c2ind,:]
    c̃[(m₁+1):end]    = qp.c₂[c2ind]

    new{T,M}(qp.Q, qp.q₁, qp.q₂, qp.A, qp.b, C̃, c̃, qp.lb, qp.ub, qp.vars, qp.name)
  end

  function (::Type{CanonicalDescription})(qp::MPSDescription{T,M}) where {T<:AbstractFloat,M<:AbstractSparseMatrix}
    c1ind   = isfinite.(qp.c₁)
    c2ind   = isfinite.(qp.c₂)
    m₁, m₂  = sum(c1ind), sum(c2ind)
    m       = m₁ + m₂
    n       = length(qp.vars)
    C̃ᵢ      = Vector{Int}()
    C̃ⱼ      = Vector{Int}()
    C̃ᵥ      = Vector{T}()
    c̃       = zeros(T, m)

    rows = rowvals(qp.C)
    vals = nonzeros(qp.C)
    for col = 1:size(qp.C,2)
      for j in nzrange(qp.C, col)
        row = rows[j]
        if c1ind[row]
          idx = count(i->i,c1ind[1:row])
          push!(C̃ᵢ,idx)
          push!(C̃ⱼ,col)
          push!(C̃ᵥ,-vals[j])
          c̃[idx] = -qp.c₁[row]
        end
        if c2ind[row]
          idx = count(i->i,c2ind[1:row])
          push!(C̃ᵢ,m₁+idx)
          push!(C̃ⱼ,col)
          push!(C̃ᵥ,vals[j])
          c̃[m₁+idx] = qp.c₂[row]
        end
      end
    end
    C̃ = sparse(C̃ᵢ, C̃ⱼ, C̃ᵥ, m, n)

    new{T,M}(qp.Q, qp.q₁, qp.q₂, qp.A, qp.b, C̃, c̃, qp.lb, qp.ub, qp.vars, qp.name)
  end
end

"""
    canonical(mps::MPSDescription) -> canon

Transform the `mps` object as returned from `parseqps` into `canon`, an object
of type `CanonicalDescription`.

`canon` contains the description of the following optimization problem formulation:

```
minimize    0.5*x'*Q*x + q₁'*x + q₂
subject to       A*x = b
                 C̃*x ≤ c̃
            lb ≤  x  ≤ ub
```

where all the matrices and vectors are of type `T`, where `T` is the first
argument in `parseqps`.

See also: `name`, `variables`, `objective`, `equalities`, `inequalities`,
  `bounds` and `parseqps`.
"""
canonical(qp::MPSDescription{T}) where {T<:AbstractFloat}     =
  CanonicalDescription(qp)

name(qp::CanonicalDescription)          = qp.name
variables(qp::CanonicalDescription)     = qp.vars
objective(qp::CanonicalDescription)     = (qp.Q, qp.q₁, qp.q₂ )
equalities(qp::CanonicalDescription)    = (qp.A, qp.b         )
inequalities(qp::CanonicalDescription)  = (qp.C̃, qp.c̃         )
bounds(qp::CanonicalDescription)        = (qp.lb, qp.ub       )

function _compact(stream, ::MIME"text/plain", qp::CanonicalDescription)
  print(stream, "QP")
end

function _full(stream, ::MIME"text/plain", qp::CanonicalDescription)
  println(stream, _header(qp))
  println(stream, "---")
  println(stream, "  minimize    0.5*x'*Q*x + q₁'*x + q₂")
  println(stream, "  subject to       A*x = b")
  println(stream, "                   C̃*x ≤ c̃")
  println(stream, "              lb ≤  x  ≤ ub")
  print(stream, "---")
end

_header(qp::CanonicalDescription) = string("QP instance: ", name(qp), ".")
_header(qp::CanonicalDescription{T,M}) where {T<:AbstractFloat,M<:AbstractSparseMatrix} = string("Sparse QP instance: ", name(qp), ".")

show(stream::IO, qp::CanonicalDescription)                          =
  show(stream, MIME("text/plain"), qp)
show(stream::IO, mime::MIME"text/plain", qp::CanonicalDescription)  =
  get(stream, :compact, false) ? _compact(stream, mime, qp) : _full(stream, mime, qp)

immutable CanonicalDescription{T<:AbstractFloat}
  Q::Matrix{T}
  q₁::Vector{T}
  q₂::T
  A::Matrix{T}
  b::Vector{T}
  C̃::Matrix{T}
  c̃::Vector{T}
  lb::Vector{T}
  ub::Vector{T}
  vars::Vector{Symbol}
  name::String

  function (::Type{CanonicalDescription}){T<:AbstractFloat}(qp::MPSDescription{T})
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

    new{T}(qp.Q, qp.q₁, qp.q₂, qp.A, qp.b, C̃, c̃, qp.lb, qp.ub, qp.vars, qp.name)
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
canonical{T<:AbstractFloat}(qp::MPSDescription{T})            =
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
  println(stream, "QP instance: ", name(qp), ".")
  println(stream, "---")
  println(stream, "  minimize    0.5*x'*Q*x + q₁'*x + q₂")
  println(stream, "  subject to       A*x = b")
  println(stream, "                   C̃*x ≤ c̃")
  println(stream, "              lb ≤  x  ≤ ub")
  print(stream, "---")
end

show(stream::IO, qp::CanonicalDescription)                          =
  show(stream, MIME("text/plain"), qp)
show(stream::IO, mime::MIME"text/plain", qp::CanonicalDescription)  =
  get(stream, :compact, false) ? _compact(stream, mime, qp) : _full(stream, mime, qp)

struct CanonicalDescription{T<:AbstractFloat}
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

  function CanonicalDescription{T}(qp::MPSDescription{T}) where T<:AbstractFloat
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

CanonicalDescription(qp::MPSDescription{T}) where T<:AbstractFloat =
  CanonicalDescription{T}(qp)
canonical(qp::MPSDescription{T})            where T<:AbstractFloat =
  CanonicalDescription{T}(qp)

objective(qp::CanonicalDescription)     = (qp.Q, qp.q₁, qp.q₂ )
inequalities(qp::CanonicalDescription)  = (qp.C̃, qp.c̃         )
equalities(qp::CanonicalDescription)    = (qp.A, qp.b         )
bounds(qp::CanonicalDescription)        = (qp.lb, qp.ub       )
variables(qp::CanonicalDescription)     = qp.vars
name(qp::CanonicalDescription)          = qp.name

struct QPInstance{T<:AbstractFloat}
  Q::Matrix{T}
  q₁::Vector{T}
  q₂::T
  c₁::Vector{T}
  C::Matrix{T}
  c₂::Vector{T}
  A::Matrix{T}
  b::Vector{T}
  x₁::Vector{T}
  x₂::Vector{T}
  vars::Vector{Symbol}
  name::String
end

function QPInstance(::Type{T}, n::Int, m₁::Int, m₂::Int,
  name::AbstractString = "QP") where T<:AbstractFloat
  QPInstance{T}(zeros(T, n, n), zeros(T, n), zero(T), fill(convert(T, -Inf), m₁),
    zeros(T, m₁, n), fill(convert(T, Inf), m₁), zeros(T, m₂, n), zeros(T, m₂),
    zeros(T, n), fill(convert(T, Inf), n), Symbol[Symbol(:x_,k) for k in 1:n],
    name)
end

QPInstance(n::Int, m₁::Int, m₂::Int, name::AbstractString = "QP")           =
  QPInstance(Float64, n, m₁, m₂, name)
QPInstance(t::Type{T}, name::AbstractString = "QP") where T<:AbstractFloat  =
  QPInstance(t, 0, 0, 0, name)
QPInstance(name::AbstractString = "QP")                                     =
  QPInstance(Float64, name)

function QPInstance(::Type{T}, Q::AbstractMatrix, q₁::AbstractVector, q₂::Real,
  c₁::AbstractVector, C::AbstractMatrix, c₂::AbstractVector, A::AbstractMatrix,
  b::AbstractVector, x₁::AbstractVector, x₂::AbstractVector,
  vars::AbstractVector{S}, name::AbstractString) where T<:AbstractFloat where
  S<:Union{Symbol,Char,AbstractString}

  n₁, n₂  = size(Q)
  m₁, n₃  = size(C)
  m₂, n₄  = size(A)
  n₅      = length(vars)

  if n₁ ≠ n₂
    warn("QPInstance: `Q` matrix must be symmetric")
    throw(DomainError())
  elseif n₁ ≠ length(q₁)
    warn("QPInstance: `q₁` must have $(n₁) elements")
    throw(DomainError())
  elseif n₁ ≠ n₃
    warn("QPInstance: `C` must have $(n₁) columns")
    throw(DomainError())
  elseif m₁ ≠ length(c₁) || m₁ ≠ length(c₂)
    warn("QPInstance: Both `c₁` and `c₂` must have $(m₁) elements")
    throw(DomainError())
  elseif n₁ ≠ n₄
    warn("QPInstance: `A` must have $(n₁) columns")
    throw(DomainError())
  elseif m₂ ≠ length(b)
    warn("QPInstance: `b` must have $(m₂) elements")
    throw(DomainError())
  elseif n₁ ≠ length(x₁) || n₁ ≠ length(x₂)
    warn("QPInstance: Both `x₁` and `x₂` must have $(n₁) elements")
    throw(DomainError())
  elseif n₁ ≠ length(vars)
    warn("QPInstance: `vars` must have $(n₁) elements")
    throw(DomainError())
  end

  QPInstance{T}(Q, q₁, q₂, c₁, C, c₂, A, b, x₁, x₂, vars, name)

end

QPInstance(Q::AbstractMatrix, q₁::AbstractVector, q₂::Real, c₁::AbstractVector,
  C::AbstractMatrix, c₂::AbstractVector, A::AbstractMatrix, b::AbstractVector,
  x₁::AbstractVector, x₂::AbstractVector, vars::AbstractVector{S},
  name::AbstractString) where S<:Union{Symbol,Char,AbstractString} =
  QPInstance(Float64, Q, q₁, q₂, c₁, C, c₂, A, b, x₁, x₂, vars, name)


objective(qp::QPInstance)     = (qp.Q, qp.q₁, qp.q₂)
inequalities(qp::QPInstance)  = (qp.C, qp.c₁, qp.c₂)
equalities(qp::QPInstance)    = (qp.A, qp.b)
bounds(qp::QPInstance)        = (qp.x₁, qp.x₂)
variables(qp::QPInstance)     = qp.vars
name(qp::QPInstance)          = qp.name

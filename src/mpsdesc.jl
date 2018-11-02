struct MPSDescription{T<:AbstractFloat,M<:AbstractMatrix}
  Q::M
  q₁::Vector{T}
  q₂::T
  c₁::Vector{T}
  C::M
  c₂::Vector{T}
  A::M
  b::Vector{T}
  lb::Vector{T}
  ub::Vector{T}
  vars::Vector{Symbol}
  name::String
end

function (::Type{MPSDescription})(::Type{T}, ::Type{Matrix}, n::Int, m₁::Int,
  m₂::Int, name::AbstractString = "QP") where {T<:AbstractFloat}
  MPSDescription{T,Matrix{T}}(zeros(T, n, n), zeros(T, n), zero(T), fill(convert(T, -Inf), m₁),
    zeros(T, m₁, n), fill(convert(T, Inf), m₁), zeros(T, m₂, n), zeros(T, m₂),
    zeros(T, n), fill(convert(T, Inf), n), Symbol[Symbol(:x_,k) for k in 1:n],
    name)
end

function (::Type{MPSDescription})(::Type{T}, ::Type{SparseMatrixCSC}, n::Int, m₁::Int,
  m₂::Int, name::AbstractString = "QP") where {T<:AbstractFloat}
  MPSDescription{T,SparseMatrixCSC{T,Int}}(spzeros(T, n, n), zeros(T, n), zero(T), fill(convert(T, -Inf), m₁),
    spzeros(T, m₁, n), fill(convert(T, Inf), m₁), spzeros(T, m₂, n), zeros(T, m₂),
    zeros(T, n), fill(convert(T, Inf), n), Symbol[Symbol(:x_,k) for k in 1:n],
    name)
end

MPSDescription(n::Int, m₁::Int, m₂::Int, name::AbstractString = "QP")     =
  MPSDescription(Float64, Matrix, n, m₁, m₂, name)
MPSDescription(t::Type{T}, name::AbstractString = "QP") where {T<:AbstractFloat} =
  MPSDescription(t, 0, 0, 0, name)
MPSDescription(name::AbstractString = "QP")                               =
  MPSDescription(Float64, name)

SparseMPSDescription(n::Int, m₁::Int, m₂::Int, name::AbstractString = "QP")     =
  MPSDescription(Float64, SparseMatrixCSC, n, m₁, m₂, name)
SparseMPSDescription(t::Type{T}, name::AbstractString = "QP") where {T<:AbstractFloat} =
  SparseMPSDescription(t, 0, 0, 0, name)
SparseMPSDescription(name::AbstractString = "QP")                               =
  SparseMPSDescription(Float64, name)

function MPSDescription(::Type{T}, Q::AbstractMatrix, q₁::AbstractVector, q₂::Real,
  c₁::AbstractVector, C::AbstractMatrix, c₂::AbstractVector, A::AbstractMatrix, b::AbstractVector,
  lb::AbstractVector, ub::AbstractVector, vars::AbstractVector{S},
  name::AbstractString) where {T<:AbstractFloat,S<:Union{Symbol,Char,AbstractString}}

  Q, C, A = convert.(Matrix{T}, (Q, C, A))

  _mpscheck(Q, q₁, q₂, c₁, C, c₂, A, b, lb, ub, vars)

  MPSDescription{T,Matrix{T}}(Q, q₁, q₂, c₁, C, c₂, A, b, lb, ub, Symbol.(vars), name)
end

function MPSDescription(::Type{T}, Q::AbstractSparseMatrix, q₁::AbstractVector, q₂::Real,
  c₁::AbstractVector, C::AbstractSparseMatrix, c₂::AbstractVector, A::AbstractSparseMatrix, b::AbstractVector,
  lb::AbstractVector, ub::AbstractVector, vars::AbstractVector{S},
  name::AbstractString) where {T<:AbstractFloat,S<:Union{Symbol,Char,AbstractString}}

  Q, C, A = convert.(SparseMatrixCSC{T,Int}, (Q, C, A))

  _mpscheck(Q, q₁, q₂, c₁, C, c₂, A, b, lb, ub, vars)

  MPSDescription{T,SparseMatrixCSC{T,Int}}(Q, q₁, q₂, c₁, C, c₂, A, b, lb, ub, Symbol.(vars), name)
end

function _mpscheck(Q::AbstractMatrix, q₁::AbstractVector, q₂::Real, c₁::AbstractVector,
  C::AbstractMatrix, c₂::AbstractVector, A::AbstractMatrix, b::AbstractVector,
  lb::AbstractVector, ub::AbstractVector, vars::AbstractVector{S}) where {S<:Union{Symbol,Char,AbstractString}}

  n₁, n₂  = size(Q)
  m₁, n₃  = size(C)
  m₂, n₄  = size(A)
  n₅      = length(vars)

  if n₁ ≠ n₂
    throw(DomainError(Q ,"MPSDescription: `Q` matrix must be symmetric"))
  elseif n₁ ≠ length(q₁)
    throw(DomainError(q₁ ,"MPSDescription: `q₁` must have $(n₁) elements"))
  elseif n₁ ≠ n₃
    throw(DomainError(C ,"MPSDescription: `C` must have $(n₁) columns"))
  elseif m₁ ≠ length(c₁) || m₁ ≠ length(c₂)
    throw(DomainError((c₁, c₂) ,"MPSDescription: Both `c₁` and `c₂` must have $(m₁) elements"))
  elseif n₁ ≠ n₄
    throw(DomainError(A ,"MPSDescription: `A` must have $(n₁) columns"))
  elseif m₂ ≠ length(b)
    throw(DomainError(b ,"MPSDescription: `b` must have $(m₂) elements"))
  elseif n₁ ≠ length(lb) || n₁ ≠ length(ub)
    throw(DomainError((lb, ub) ,"MPSDescription: Both `lb` and `ub` must have $(n₁) elements"))
  elseif n₁ ≠ length(vars)
    throw(DomainError(vars ,"MPSDescription: `vars` must have $(n₁) elements"))
  end
end

MPSDescription(Q::AbstractMatrix, q₁::AbstractVector, q₂::Real,
  c₁::AbstractVector, C::AbstractMatrix, c₂::AbstractVector,
  A::AbstractMatrix, b::AbstractVector, lb::AbstractVector,
  ub::AbstractVector, vars::AbstractVector{S}, name::AbstractString) where {S<:Union{Symbol,Char,AbstractString}} =
  MPSDescription(Float64, Q, q₁, q₂, c₁, C, c₂, A, b, lb, ub, vars, name)

"""
    name(qp) -> qpname

Return the problem instance `qp`'s name. `qp` can be either `MPSDescription` or
`CanonicalDescription`. `qpname` is a `String`.

See also: `variables`, `objective`, `equalities`, `inequalities`, `bounds`,
  `parseqps` and `canonical`.
"""
name(qp::MPSDescription)          = qp.name

"""
    variables(qp) -> vars

Return the problem instance `qp`'s variable names. `qp` can be either
`MPSDescription` or `CanonicalDescription`. `vars` is a `Vector{Symbol}`.

See also: `name`, `objective`, `equalities`, `inequalities`, `bounds`, `parseqps`
  and `canonical`.
"""
variables(qp::MPSDescription)     = qp.vars

"""
    objective(qp) -> (Q, q₁, q₂)

Return the problem instance `qp`'s objective part. `qp` can be either
`MPSDescription` or `CanonicalDescription`. `Q` is a `Matrix{T}`, `q₁` is a
`Vector{T}` and `q₂` is a `T`, where `T` is the first argument in `parseqps`.

See also: `name`, `variables`, `equalities`, `inequalities`, `bounds`, `parseqps`
  and `canonical`.
"""
objective(qp::MPSDescription)     = (qp.Q, qp.q₁, qp.q₂ )

"""
    equalities(qp) -> (A, b)

Return the problem instance `qp`'s `A, b` pair that defines the equalities among
the problem's variables. `qp` can be either `MPSDescription` or `CanonicalDescription`.
`A` is a `Matrix{T}` and `b` is a `Vector{T}`, where `T` is the first argument
in `parseqps`.

See also: `name`, `variables`, `objective`, `inequalities`, `bounds`, `parseqps`
  and `canonical`.
"""
equalities(qp::MPSDescription)    = (qp.A, qp.b         )

"""
    inequalities(qp::MPSDescription)        -> (C, c₁, c₂)
    inequalities(qp::CanonicalDescription)  -> (C̃, c̃)

Return the problem instance `qp`'s part that defines the inequalities among
the problem's variables.

When `qp` is `MPSDescription`, `C` is a `Matrix{T}`, and `c₁` and `c₂` are
`Vector{T}`, where `T` is the first argument in `parseqps`. When `qp` is
`CanonicalDescription`, `C̃` is a `Matrix{T}` and `c̃` is a `Vector{T}`.

See also: `name`, `variables`, `objective`, `equalities`, `bounds`, `parseqps`
and `canonical`.
"""
inequalities(qp::MPSDescription)  = (qp.C, qp.c₁, qp.c₂ )

"""
    bounds(qp) -> (lb, ub)

Return the lower and upper bounds, `lb` and `ub`, defined on the variables of the
problem instance `qp`. `qp` can be either `MPSDescription` or `CanonicalDescription`.
Both `lb` and `ub` are `Vector{T}`, where `T` is the first argument in `parseqps`.

See also: `name`, `variables`, `objective`, `inequalities`, `inequalities`,
  `parseqps` and `canonical`.
"""
bounds(qp::MPSDescription)        = (qp.lb, qp.ub       )

function _compact(stream, ::MIME"text/plain", qp::MPSDescription)
  print(stream, "QP")
end

function _full(stream, ::MIME"text/plain", qp::MPSDescription)
  println(stream, _header(qp))
  println(stream, "---")
  println(stream, "  minimize    0.5*x'*Q*x + q₁'*x + q₂")
  println(stream, "  subject to       A*x = b")
  println(stream, "              c₁ ≤ C*x ≤ c₂")
  println(stream, "              lb ≤  x  ≤ ub")
  print(stream, "---")
end

_header(qp::MPSDescription) = string("QP instance: ", name(qp), ".")
_header(qp::MPSDescription{T,M}) where {T<:AbstractFloat,M<:AbstractSparseMatrix} =
  string("Sparse QP instance: ", name(qp), ".")

show(stream::IO, qp::MPSDescription)                          =
  show(stream, MIME("text/plain"), qp)
show(stream::IO, mime::MIME"text/plain", qp::MPSDescription)  =
  get(stream, :compact, false) ? _compact(stream, mime, qp) : _full(stream, mime, qp)

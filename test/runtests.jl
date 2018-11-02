using QPSParser
using QPSParser: MPSDescription
using Test, LinearAlgebra, SparseArrays

import InteractiveUtils: subtypes

eye(n,m) = Matrix{Float64}(I, n, m)
eye(n) = eye(n,n)

speye(n,m) = SparseMatrixCSC{Float64,Int}(I,n,m)
speye(n) = speye(n,n)

include("mpsdesc.jl")
include("parseqps.jl")

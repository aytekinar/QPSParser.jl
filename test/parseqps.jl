@testset "QP Example from Maros" begin
  @testset "parseqps($T, filename)" for T in subtypes(AbstractFloat)
    # MPS description
    mps   = parseqps(T, joinpath(dirname(@__FILE__), "qp-example.qps"))

    @test name(mps) == "QP example"

    @test variables(mps) == [:c1, :c2]

    Q, q₁, q₂ = objective(mps)
    @test Q ≈ [8 2; 2 10] && eltype(Q) == T
    @test q₁ ≈ [1.5; -2] && eltype(q₁) == T
    @test q₂ ≈ 4 && eltype(q₂) == T

    A, b = equalities(mps)
    @test isempty(A) == true && eltype(A) == T
    @test isempty(b) == true && eltype(b) == T

    C, c₁, c₂ = inequalities(mps)
    @test C ≈ [2 1; -1 2] && eltype(C) == T
    @test c₁ ≈ [2; -Inf] && eltype(c₁) == T
    @test c₂ ≈ [Inf, 6] && eltype(c₂) == T

    lb, ub = bounds(mps)
    @test lb ≈ zeros(2) && eltype(lb) == T
    @test ub ≈ [20, Inf] && eltype(ub) == T
    # Canonical description
    canon = canonical(mps)

    @test name(canon) == name(mps)

    @test variables(canon) == variables(mps)

    QQ, qq1, qq2 = objective(canon)
    @test QQ ≈ Q && qq1 ≈ q₁ && qq2 ≈ q₂

    AA, bb = equalities(canon)
    @test AA ≈ A && bb ≈ b

    C̃, c̃ = inequalities(canon)
    @test C̃ ≈ [-2 -1; -1 2] && eltype(C̃) == T
    @test c̃ ≈ [-2; 6] && eltype(c̃) == T

    llb, uub = bounds(canon)
    @test llb ≈ lb && uub ≈ ub
  end

  @testset "parse_sparseqps($T, filename)" for T in subtypes(AbstractFloat)
    # MPS description
    mps   = parse_sparseqps(T, joinpath(dirname(@__FILE__), "qp-example.qps"))

    @test name(mps) == "QP example"

    @test variables(mps) == [:c1, :c2]

    Q, q₁, q₂ = objective(mps)
    @test typeof(Q) == SparseMatrixCSC{T,Int} && Matrix(Q) ≈ [8 2; 2 10]
    @test q₁ ≈ [1.5; -2] && eltype(q₁) == T
    @test q₂ ≈ 4 && eltype(q₂) == T

    A, b = equalities(mps)
    @test typeof(A) == SparseMatrixCSC{T,Int} && isempty(A) == true
    @test isempty(b) == true && eltype(b) == T

    C, c₁, c₂ = inequalities(mps)
    @test typeof(C) == SparseMatrixCSC{T,Int} && Matrix(C) ≈ [2 1; -1 2]
    @test c₁ ≈ [2; -Inf] && eltype(c₁) == T
    @test c₂ ≈ [Inf, 6] && eltype(c₂) == T

    lb, ub = bounds(mps)
    @test lb ≈ zeros(2) && eltype(lb) == T
    @test ub ≈ [20, Inf] && eltype(ub) == T
    # Canonical description
    canon = canonical(mps)

    @test name(canon) == name(mps)

    @test variables(canon) == variables(mps)

    QQ, qq1, qq2 = objective(canon)
    @test QQ ≈ Q && qq1 ≈ q₁ && qq2 ≈ q₂

    AA, bb = equalities(canon)
    @test AA ≈ A && bb ≈ b

    C̃, c̃ = inequalities(canon)
    @test typeof(C̃) == SparseMatrixCSC{T,Int} && Matrix(C̃) ≈ [-2 -1; -1 2]
    @test c̃ ≈ [-2; 6] && eltype(c̃) == T

    llb, uub = bounds(canon)
    @test llb ≈ lb && uub ≈ ub
  end
end

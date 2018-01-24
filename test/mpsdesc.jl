@testset "MPS Description" for (Q, C, A) in [(eye(2),eye(3,2),eye(4,2)),(speye(2),speye(3,2),speye(4,2))]
  q₁      = ones(2)
  q₂      = 0.
  c₁      = zeros(3)
  c₂      = ones(3)
  b       = ones(4)
  lb      = zeros(2)
  ub      = ones(2)
  vars    = ["x₁", "x₂"]
  qpname  = "Test Problem"
  @testset "Domain Errors" begin
    @test_throws DomainError MPSDescription(eye(3,2), q₁, q₂, c₁, C, c₂, A, b, lb, ub, vars, qpname)
    @test_throws DomainError MPSDescription(Q, ones(3), q₂, c₁, C, c₂, A, b, lb, ub, vars, qpname)
    @test_throws DomainError MPSDescription(Q, q₁, q₂, c₁, eye(3), c₂, A, b, lb, ub, vars, qpname)
    @test_throws DomainError MPSDescription(Q, q₁, q₂, zeros(2), C, c₂, A, b, lb, ub, vars, qpname)
    @test_throws DomainError MPSDescription(Q, q₁, q₂, c₁, C, ones(2), A, b, lb, ub, vars, qpname)
    @test_throws DomainError MPSDescription(Q, q₁, q₂, c₁, C, c₂, eye(4,3), b, lb, ub, vars, qpname)
    @test_throws DomainError MPSDescription(Q, q₁, q₂, c₁, C, c₂, A, ones(3), lb, ub, vars, qpname)
    @test_throws DomainError MPSDescription(Q, q₁, q₂, c₁, C, c₂, A, b, zeros(3), ub, vars, qpname)
    @test_throws DomainError MPSDescription(Q, q₁, q₂, c₁, C, c₂, A, b, lb, ones(3), vars, qpname)
    @test_throws DomainError MPSDescription(Q, q₁, q₂, c₁, C, c₂, A, b, lb, ub, ["x₁", "x₂", "x₃"], qpname)
  end
  @testset "Methods for MPSDescription{$T}" for T in subtypes(AbstractFloat)
    problem = MPSDescription(T, Q, q₁, q₂, c₁, C, c₂, A, b, lb, ub, vars, qpname)

    @test name(problem) == qpname

    @test variables(problem) == Symbol.(vars)

    QQ, qq1, qq2 = objective(problem)
    @test QQ ≈ Q && eltype(QQ) == T
    @test qq1 ≈ q₁ && eltype(qq1) == T
    @test qq2 ≈ q₂ && eltype(qq2) == T

    AA, bb = equalities(problem)
    @test AA ≈ A && eltype(AA) == T
    @test bb ≈ b && eltype(bb) == T

    CC, cc1, cc2 = inequalities(problem)
    @test CC ≈ C && eltype(CC) == T
    @test cc1 ≈ c₁ && eltype(cc1) == T
    @test cc2 ≈ c₂ && eltype(cc2) == T

    llb, uub = bounds(problem)
    @test llb ≈ lb && eltype(llb) == T
    @test uub ≈ ub && eltype(uub) == T
  end
end

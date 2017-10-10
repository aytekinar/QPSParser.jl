# QPSParser

[![Unix][unix-img]][unix-link] [![Windows][win-img]][win-link]
[![Coveralls][ca-img]][ca-link] [![Codecov][cc-img]][cc-link]

Julia package for parsing QPS files.

[unix-img]: https://img.shields.io/travis/aytekinar/QPSParser.jl/master.svg?label=unix
[unix-link]: https://travis-ci.org/aytekinar/QPSParser.jl
[win-img]: https://img.shields.io/appveyor/ci/aytekinar/qpsparser-jl/master.svg?label=windows
[win-link]: https://ci.appveyor.com/project/aytekinar/qpsparser-jl/branch/master
[ca-img]: https://img.shields.io/coveralls/aytekinar/QPSParser.jl/master.svg?label=coveralls
[ca-link]: https://coveralls.io/github/aytekinar/QPSParser.jl?branch=master
[cc-img]: https://img.shields.io/codecov/c/github/aytekinar/QPSParser.jl/master.svg?label=codecov
[cc-link]: https://codecov.io/gh/aytekinar/QPSParser.jl?branch=master

## Description

QPS is a file format, which is proposed by [Maros and Mészáros][maros-doi], used
to encode convex quadratic programming problems in plain text files. It is an
extension to the well-known [MPS format][mps-link]. The authors also provide a
set of large and/or ill-conditioned QP test problems on their [website][test-set]
to be used to assess the performance of different solvers.

`QPSParser` is a Julia package that aims at helping users parse problems written
in the QPS file format.

[maros-doi]: http://dx.doi.org/10.1080/10556789908805768
[mps-link]: https://en.wikipedia.org/wiki/MPS_(format)
[test-set]: http://www.doc.ic.ac.uk/~im/

## Basic Usage

After installing the package via `Pkg.add("QPSParser")`, we import the useful
parts of the package to our workspace by `using QPSParser`.

### Parsing QPS files

To parse a QPS file, we use `parseqps("filename.qps")`. The output will be an
object of type `MPSDescription`, which is the description of the following
optimization problem as defined in `filename.qps`:

```
minimize    0.5*x'*Q*x + q₁'*x + q₂
subject to       A*x = b
            c₁ ≤ C*x ≤ c₂
            lb ≤  x  ≤ ub
```

To illustrate, let's assume that the `QP Example` given in [Maros and Mészáros][maros-doi]

```
NAME          QP example
ROWS
 N  obj
 G  r1
 L  r2
COLUMNS
    c1        r1                2.0    r2               -1.0
    c1        obj               1.5
    c2        r1                1.0    r2                2.0
    c2        obj              -2.0
RHS
    rhs1      obj              -4.0
    rhs1      r1                2.0    r2                6.0
RANGES
BOUNDS
UP  bnd1      c1               20.0
QUADOBJ
    c1        c1                8.0
    c1        c2                2.0
    c2        c2               10.0
ENDATA
```

is saved as `filename.qps`. Then,

```julia
julia> mps   = parseqps("test/QPExample.qps")
QP instance: QP example.
---
  minimize    0.5*x'*Q*x + q₁'*x + q₂
  subject to       A*x = b
              c₁ ≤ C*x ≤ c₂
              lb ≤  x  ≤ ub
---
```

For more information on the function, type `?parseqps` in Julia.

### Accessing the Fields

There are convenient getter functions defined for the `MPSDescription` objects.

For instance, to get the `name` and `variables` of the given problem, we can
write

```julia
julia> name(mps)
"QP example"

julia> variables(mps)
2-element Array{Symbol,1}:
 :c1
 :c2
```

Then, we can get the `objective` part of the example problem description above
by typing

```julia
julia> Q, q₁, q₂ = objective(mps);

julia> Q
2×2 Array{Float64,2}:
 8.0   2.0
 2.0  10.0

julia> q₁
2-element Array{Float64,1}:
  1.5
 -2.0

julia> q₂
4.0
```

Similarly, to get the `equalities` and `inequalities` defined on the `variables`
of the problem, we can write

```julia
julia> A, b = equalities(mps);

julia> A  # no equalities given in the example
0×2 Array{Float64,2}

julia> b  # no equalities given in the example
0-element Array{Float64,1}

julia> C, c₁, c₂ = inequalities(mps);

julia> C
2×2 Array{Float64,2}:
  2.0  1.0
 -1.0  2.0

julia> c₁
2-element Array{Float64,1}:
    2.0
 -Inf

julia> c₂
2-element Array{Float64,1}:
 Inf
   6.0
```

Finally, the lower and upper `bounds` defined on the `variables` of the problem
can be obtained by

```julia
julia> lb, ub = bounds(mps);

julia> lb
2-element Array{Float64,1}:
 0.0
 0.0

julia> ub
2-element Array{Float64,1}:
  20.0
 Inf
```

Please note that, in `MPSDescription` objects, `c₁`, `c₂`, `lb` and `ub` can all
contain infinite values.

### Obtaining a Canonical Representation

Sometimes, it might be useful to obtain the following canonical representation
of the same optimization problem defined above:

```
minimize    0.5*x'*Q*x + q₁'*x + q₂
subject to       A*x = b
                 C̃*x ≤ c̃
            lb ≤  x  ≤ ub
```

To achieve this, we can write `canon = canonical(mps)`, which returns a
`CanonicalDescription` object:

```julia
julia> canon = canonical(mps)
QP instance: QP example.
---
  minimize    0.5*x'*Q*x + q₁'*x + q₂
  subject to       A*x = b
                   C̃*x ≤ c̃
              lb ≤  x  ≤ ub
---

julia> name(canon);                   # same as `name(mps)`

julia> variables(canon);              # same as `variables(mps)`

julia> Q, q₁, q₂ = objective(canon);  # same as `objective(mps)`

julia> A, b = equalities(canon);      # same as `equalities(mps)`

julia> lb, ub = bounds(canon);        # same as `bounds(mps)`

julia> C̃, c̃ = inequalities(canon);

julia> C̃
2×2 Array{Float64,2}:
 -2.0  -1.0
 -1.0   2.0

julia> c̃
2-element Array{Float64,1}:
 -2.0
  6.0
```

## Limitations

`QPSParser`,

- does **not** guarantee that `Q` in the problem descriptions is positive
  semi-definite. Any user of the package *has to* check the definiteness of
  `Q` in their applications,
- gives both `MPSDescription` and `CanonicalDescription` in **dense** format.
  In other words, all the matrices and vectors involved in the descriptions
  are dense matrices and vectors. This might be a limitation when dealing with
  large problems or using the package in low-memory platforms.

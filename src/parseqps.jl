"""
    parseqps([t::Type{T}, ]filename::AbstractString) -> mps

Parse QPS formatted file `filename`. `T` must be an `AbstractFloat` type. When
only `filename` is provided, the first argument defaults to `Float64`.

`mps` is of `MPSDescription` type. It contains the description of the the
following optimization problem formuliation:

```
minimize    0.5*x'*Q*x + q₁'*x + q₂
subject to       A*x = b
            c₁ ≤ C*x ≤ c₂
            lb ≤  x  ≤ ub
```

where all the matrices and vectors are of type `T`.

See also: `name`, `variables`, `objective`, `equalities`, `inequalities`,
  `bounds` and `canonical`.
"""
function parseqps{T<:AbstractFloat}(::Type{T}, filename::AbstractString)
  # Initialize auxiliary variables
  qpname                        = "QP"
  mode                          = :NAME
  ineqrowidx, eqrowidx, varidx  = 0, 0, 1
  objgiven                      = false
  objsymbol                     = :obj

  # Define sections
  vars    = Dict{Symbol,Int}()
  rows    = Dict{Symbol,Pair{Int,Symbol}}()
  cols    = Dict{Symbol,Vector{Pair{Symbol,T}}}()
  rhs     = Dict{Symbol,T}()
  ranges  = Dict{Symbol,T}()
  bounds  = Dict{Symbol,Vector{Pair{Symbol,T}}}()
  quadobj = Dict{Symbol,Vector{Pair{Symbol,T}}}()

  # Parse the file
  open(filename) do io
    firstline = split(readline(io))
    if Symbol(firstline[1]) == :NAME
      qpname = join(firstline[2:end], " ")
    else
      throw(ArgumentError("`NAME` field is expected on the first line"))
    end
    for line in eachline(io)
      words = split(line)
      length(words) == 1 && (mode = Symbol(words[1]); continue)
      if mode == :ROWS
        rowsym = Symbol(words[2])
        rows[rowsym] =
          words[1] == "N" ? (objgiven = true; objsymbol = rowsym; Pair(0, :obj)):
          words[1] == "L" ? (ineqrowidx += 1; Pair(ineqrowidx, :leq))           :
          words[1] == "G" ? (ineqrowidx += 1; Pair(ineqrowidx, :geq))           :
                            (eqrowidx   += 1; Pair(eqrowidx, :eq))
      elseif mode == :COLUMNS
        var = Symbol(words[1])
        if get!(vars, var, 0) == 0
          vars[var] = varidx
          varidx    += 1
        end
        for idx = 2:2:length(words)
          push!(get!(cols, var, Pair{Symbol,T}[]),
            Pair(Symbol(words[idx]), convert(T, parse(Float64, words[idx+1]))))
        end
      elseif mode == :RHS
        # NOTE: Discarding rhs name, *i.e.*, assuming only one rhs.
        for idx = 2:2:length(words)
          rhs[Symbol(words[idx])] = convert(T, parse(Float64, words[idx+1]))
        end
      elseif mode == :RANGES
        # NOTE: Discarding range name, *i.e.*, assuming only one range.
        for idx = 2:2:length(words)
          ranges[Symbol(words[idx])] = convert(T, parse(Float64, words[idx+1]))
        end
      elseif mode == :BOUNDS
        var = Symbol(words[3])
        bnd = words[1] == "LO" ? :lower :
              words[1] == "UP" ? :upper :
              words[1] == "FR" ? :free  : :fixed
        push!(get!(bounds, var, Pair{Symbol,T}[]),
          Pair(bnd, convert(T, parse(Float64, bnd == :free ? "0" : words[4]))))
      elseif mode == :QUADOBJ
        var = Symbol(words[1])
        for idx = 2:2:length(words)
          push!(get!(quadobj, var, Pair{Symbol,T}[]),
            Pair(Symbol(words[idx]), convert(T, parse(Float64, words[idx+1]))))
        end
      elseif mode == :ENDDATA
        break
      else
        throw(ArgumentError("$(mode) is not a valid word"))
      end
    end

    # Create the matrices
    n, m₁, m₂ = length(vars), ineqrowidx, eqrowidx
    Q         = zeros(T, n, n )
    q₁        = zeros(T, n)
    q₂        = zero(T)
    c₁        = fill(convert(T, -Inf), m₁)
    C         = zeros(T, m₁, n)
    c₂        = fill(convert(T, Inf), m₁)
    A         = zeros(T, m₂, n)
    b         = zeros(T, m₂)
    lb        = zeros(T, n)
    ub        = fill(convert(T, Inf), n)

    varvec    = Vector{Symbol}(n)
    for (varname, varidx) in vars
      varvec[varidx] = varname
    end

    # Fill (in)equalities with ranges
    for (key, val) in cols
      colidx = vars[key]
      for (rowsymbol, rowval) in val
        rowidx, rowtype = rows[rowsymbol]
        if rowidx == 0
          q₁[colidx] = rowval
        elseif rowtype == :eq
          A[rowidx, colidx] = rowval
          b[rowidx]         = get(rhs, rowsymbol, zero(T))
        elseif rowtype == :leq
          C[rowidx, colidx] = rowval
          c₂[rowidx]        = get(rhs, rowsymbol, zero(T))
          if haskey(ranges, rowsymbol)
            c₁[rowidx]      = c₂[rowidx] - abs(ranges[rowsymbol])
          end
        else
          C[rowidx, colidx] = rowval
          c₁[rowidx]        = get(rhs, rowsymbol, zero(T))
          if haskey(ranges, rowsymbol)
            c₂[rowidx]      = c₁[rowidx] + abs(ranges[rowsymbol])
          end
        end
      end
    end

    if objgiven
      q₂ = -get(rhs, objsymbol, zero(T))
    end

    # Fill bounds
    for (key, val) in bounds
      colidx = vars[key]
      for (bndtype, bndval) in val
        if bndtype == :lower
          lb[colidx] = bndval
        elseif bndtype == :upper
          ub[colidx] = bndval
        elseif bndtype == :fixed
          lb[colidx] = ub[colidx] = bndval
        else
          lb[colidx] = convert(T, -Inf)
        end
      end
    end

    # Fill Q
    for (colkey, val) in quadobj
      colidx = vars[colkey]
      for (rowkey, rowval) in val
        rowidx = vars[rowkey]
        Q[rowidx, colidx] = Q[colidx, rowidx] = rowval
      end
    end

    # Return the problem
    return MPSDescription(T, Q, q₁, q₂, c₁, C, c₂, A, b, lb, ub, varvec, qpname)
  end
end

parseqps(filename::AbstractString) = parseqps(Float64, filename)

module QPSParse

export  MPSDescription,
        CanonicalDescription,
        canonical,
        objective,
        inequalities,
        equalities,
        bounds,
        variables,
        name,
        qpsparse

include("mpsdesc.jl")
include("canonicaldesc.jl")
include("parse.jl")

end # module

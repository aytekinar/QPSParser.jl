module QPSParse

export  MPSDescription,
        objective,
        inequalities,
        equalities,
        bounds,
        variables,
        name,
        qpsparse

include("mpsdesc.jl")
include("parse.jl")

end # module

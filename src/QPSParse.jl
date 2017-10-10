module QPSParse

export  QPInstance,
        objective,
        inequalities,
        equalities,
        bounds,
        variables,
        name,
        qpsparse

include("QPInstance.jl")
include("parse.jl")

end # module

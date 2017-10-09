module QPS

export  QPInstance,
        objective,
        inequalities,
        equalities,
        bounds,
        variables,
        name,
        parseqps

include("QPInstance.jl")
include("parseqps.jl")

end # module

module QPSParser

import Base: show

export  canonical,
        objective,
        inequalities,
        equalities,
        bounds,
        variables,
        name,
        parseqps,
        parse_sparseqps

include("mpsdesc.jl")
include("canonicaldesc.jl")
include("parseqps.jl")

end # module

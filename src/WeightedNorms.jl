VERSION >= v"0.4.0-dev+6521" && __precompile__()
module WeightedNorms

using Compat

include("norms.jl")

include("traces.jl")

end # module

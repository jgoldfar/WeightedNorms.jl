VERSION >= v"0.4.0-dev+6521" && __precompile__()
module WeightedNorms

# Make true to emulate approximation of Lebesgue norms
const honestymode = true

using Compat

include("norms.jl")

include("traces.jl")

end # module

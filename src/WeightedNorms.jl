VERSION >= v"0.4.0-dev+6521" && __precompile__()
module WeightedNorms

@static if VERSION < v"0.7-"
    using Compat: AbstractRange
end

include("norms.jl")

include("traces.jl")

end # module

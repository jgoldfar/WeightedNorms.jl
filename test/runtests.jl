@static if VERSION < v"0.7-"
    using Base.Test
    using Compat: range
else
    using Test
end

using WeightedNorms
const WN = WeightedNorms

const verbose = false

@testset "WeightedNorms" begin
    include("norms.jl")
    include("traces.jl")
end

if VERSION < v"0.7-"
    using Base.Test
    using Compat
    const range = Compat.range
    @static if VERSION < v"0.6-"
        range{T}(start::T; stop::T=one(T), length::Int = 100) = linspace(start, stop, length)
    end
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

if VERSION >= v"0.5-"
  using Base.Test
else
  using BaseTestNext
  const Test = BaseTestNext
end
const verbose = false
const pkgbasedir = dirname(dirname(@__FILE__))
include("norms.jl")

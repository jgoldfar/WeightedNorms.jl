const testGridLength = 100
const tolerance = 1e-10
testgridArchetype = range(0, stop=1, length=testGridLength)
testGrids = [testgridArchetype, collect(testgridArchetype)]

func1(x,t) = first(promote(zero(x), zero(t)))
func2(x,t) = first(promote(one(x), one(t)))
func3(x,t) = x*t
func4(x,t) = sin(x)-cos(t)
testFuncs = [func1, func2, func3, func4]

@testset "space_trace_*, grid=$testgrid" for testgrid in testGrids
    @testset "Positivity, testFunction $func" for func in testFuncs
        U = [func(x,t) for x in testgrid, t in testgrid]
        @test WN.space_trace_norm(U, testgrid, testGridLength, testGridLength-1) >= 0

        v = [func(x, testgrid[testGridLength-1]) for x in testgrid]
        @test WN.space_trace_dist(U, v, testgrid, testGridLength, testGridLength-1) >= 0
    end
end

@testset "space_dist, grid=$testgrid" for testgrid in testGrids
    @testset "Identity, function=$func" for func in testFuncs
        U = [func(x,t) for x in testgrid, t in testgrid]
        @test isapprox(WN.space_dist(U, U, testgrid, testgrid[2] - testgrid[1]), 0)
    end
end

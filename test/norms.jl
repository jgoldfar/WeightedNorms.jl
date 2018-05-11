const testGridLength = 100
const tolerance = 1e-10
testgridArchetype = range(0, stop=1, length=testGridLength)
testGrids = [testgridArchetype, collect(testgridArchetype)]

@testset "Grid type=$(typeof(testgrid))" for testgrid in testGrids
    zerofunc(t) = zero(t)
    constfunc(t) = one(t)
    linfunc(t) = t
    sqrfunc(t) = (t^2) / 2
    cubefunc(t) = (t^3) / 6
    signchangingfunc(t) = t - 1/2
    testFunctions = [zerofunc, constfunc, linfunc, sqrfunc, cubefunc, signchangingfunc, sin, exp]

    @testset "Equality, norm=$func" for func in WN.exportedNorms
        @testset "Equality for Equally Spaced Grids, test function=$testfunc" for testfunc in testFunctions
            testdata = map(testfunc, testgrid)
            f1 = @inferred func(testdata, testgrid)
            f2 = @inferred func(testdata, testgrid[2] - testgrid[1])
            @test isapprox(f1, f2, atol = tolerance)
        end
    end

    @testset "Known Exact Norms" for func in WN.exportedNorms
        @testset "Norm of Zero is Zero" begin
            # Norm of zero is zero
            testdata = map(zerofunc, testgrid)
            f1 = @inferred func(testdata, testgrid)
            f2 = @inferred func(testdata, testgrid[2] - testgrid[1])
            @test isapprox(f1, 0.0)
            @test isapprox(f2, 0.0)
        end
        #   Norm is non-negative
        @testset "Positivity, func=$testfunc"for testfunc in testFunctions
            testdata = map(testfunc, testgrid)
            f1 = @inferred func(testdata, testgrid)
            f2 = @inferred func(testdata, testgrid[2] - testgrid[1])
            @test f1 >= 0.0
            @test f2 >= 0.0
        end
    end

    # Other known values
    knownValues = (
                   (constfunc, WN.norm_l2, 1.0),
                   (constfunc, WN.normsq_w21_only, 0.0),
                   (constfunc, WN.normsq_w22_only, 0.0),
                   (constfunc, WN.norm_w21_only, 0.0),
                   (constfunc, WN.norm_w22_only, 0.0),
                   (linfunc, WN.norm_w22_only, 0.0),
                   (linfunc, WN.normsq_w22_only, 0.0),
                   (linfunc, WN.norm_w21_only, 1.0),
                   #(sqrfunc, WN.norm_w22_only, 1.0)
                   )
    @testset "Known Values, testfunc=$testfunc, normfunc=$normfunc" for (testfunc, normfunc, expected) in knownValues
        testdata = map(testfunc, testgrid)
        f1 = @inferred normfunc(testdata, testgrid)
        f2 = @inferred normfunc(testdata, testgrid[2] - testgrid[1])
        @test (abs(f1 - expected) < tolerance)
        @test (abs(f2 - expected) < tolerance)
    end #for
    testOnesMat = ones(testGridLength, testGridLength)
    testZerosMat = zeros(testGridLength, testGridLength)
    @testset "Matrix Norm, func=$func" for func in [WN.norm_l2, WN.norm_l1]
        v1 = @inferred func(testOnesMat, testgrid, testgrid)
        @test isapprox(v1, 1.0)

        v2 = @inferred func(testZerosMat, testgrid, testgrid)
        @test isapprox(v2, 0.0)
    end
end

@testset "Edge case regression tests" begin
    twoValGrid = [0.0, 1.0]
    @test isnan(WN.norm_w22_only(twoValGrid, [-1.0, 1.0]))
    @test isnan(WN.norm_w22_only(twoValGrid, 1e-1))

    threeValGrid1 = [0.0, 1.0, 2.0]
    @test isapprox(WN.norm_w22_only(threeValGrid1, [1.0, 2.0, 3.0]), 0.0, atol=tolerance)
    @test isapprox(WN.norm_w22_only(threeValGrid1, 1), 0.0, atol=tolerance)

    threeValGrid2 = [1.0, 0.0, 1.0]
    @test WN.norm_w22_only(threeValGrid2, [1.0, 2.0, 3.0]) >= 0
    @test WN.norm_w22_only(threeValGrid2, 1) >= 0
    @test WN.norm_w22_only(-threeValGrid2, [1.0, 2.0, 3.0]) >= 0
    @test WN.norm_w22_only(-threeValGrid2, 1) >= 0

    @test WN.norm_w22_only(-threeValGrid2, 1) == WN.norm_w22_only(threeValGrid2, 1)
    @test WN.norm_w22_only(-threeValGrid2, [1.0, 2.0, 3.0]) == WN.norm_w22_only(threeValGrid2, [1.0, 2.0, 3.0])
end

@static if VERSION < v"0.7-"
    using Base.Test
    using Compat: range
else
    using Test
end

using WeightedNorms
@static if VERSION >= v"1.0-"
    if !(@isdefined WN)
        const WN = WeightedNorms
    end
else
    if !isdefined(:WN)
        const WN = WeightedNorms
    end
end

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

@testset "Convergence for Grid Input" begin
    ConvergenceNRange = [2^n for n in 3:6]
    L1ErrorValues = zeros(length(ConvergenceNRange))
    L2ErrorValues = zeros(length(ConvergenceNRange))

    @testset "Constant Function C=$C (so norm is domain size), LX=$LX, LT=$LT" for C in -1:1, LX in 3:5, LT in 1:4
        for (i,N) in enumerate(ConvergenceNRange)
            xTestGrid = collect(range(-1, stop=(-1+LX), length=N+1))
            tTestGrid = collect(range(0, stop=(0+LT), length=N))
            testfunc = C*ones(N+1, N)

            testfuncL1NormExact = abs(C)*LX*LT
            L1ErrorValues[i] = abs(testfuncL1NormExact - WN.norm_l1(testfunc, xTestGrid, tTestGrid))

            testfuncL2NormExact = sqrt(abs(C)*LX*LT)
            L2ErrorValues[i] = abs(testfuncL2NormExact - WN.norm_l2(testfunc, xTestGrid, tTestGrid))
        end

        @test all(e <= tolerance for e in L1ErrorValues)
        @test all(e <= tolerance for e in L2ErrorValues)
    end

    fill!(L1ErrorValues, 0)
    fill!(L2ErrorValues, 0)

    # Nonlinear function with changing sign.
    for (i,N) in enumerate(ConvergenceNRange)
        xTestGrid = collect(range(-1, stop=0.5, length=N+1))
        tTestGrid = collect(range(0, stop=1, length=N))
        testfunc1 = [x * (4*t - 1) for x in xTestGrid, t in tTestGrid]

        testfunc1L1NormExact = 25//32
        L1ErrorValues[i] = abs(testfunc1L1NormExact - WN.norm_l1(testfunc1, xTestGrid, tTestGrid))

        testfunc1L2NormExact = 7//8
        L2ErrorValues[i] = abs(testfunc1L2NormExact - WN.norm_l2(testfunc1, xTestGrid, tTestGrid)^2)
    end
    @test all(d <=0 for d in diff(L1ErrorValues))
    @test all(d <=0 for d in diff(L2ErrorValues))

end

@testset "Convergence for Range Input" begin
    ConvergenceNRange = [2^n for n in 3:6]
    L1ErrorValues = zeros(length(ConvergenceNRange))
    L2ErrorValues = zeros(length(ConvergenceNRange))

    @testset "Constant Function C=$C (so norm is domain size), LX=$LX, LT=$LT" for C in -1:1, LX in 3:5, LT in 1:4
        for (i,N) in enumerate(ConvergenceNRange)
            xTestGrid = range(-1, stop=(-1+LX), length=N+1)
            tTestGrid = range(0, stop=(0+LT), length=N)
            testfunc = C*ones(N+1, N)

            testfuncL1NormExact = abs(C)*LX*LT
            L1ErrorValues[i] = abs(testfuncL1NormExact - WN.norm_l1(testfunc, xTestGrid, tTestGrid))

            testfuncL2NormExact = sqrt(abs(C)*LX*LT)
            L2ErrorValues[i] = abs(testfuncL2NormExact - WN.norm_l2(testfunc, xTestGrid, tTestGrid))
        end

        @test all(e <= tolerance for e in L1ErrorValues)
        @test all(e <= tolerance for e in L2ErrorValues)
    end

    fill!(L1ErrorValues, 0)
    fill!(L2ErrorValues, 0)

    # Nonlinear function with changing sign.
    for (i,N) in enumerate(ConvergenceNRange)
        xTestGrid = range(-1, stop=0.5, length=N+1)
        tTestGrid = range(0, stop=1, length=N)
        testfunc1 = [x * (4*t - 1) for x in xTestGrid, t in tTestGrid]

        testfunc1L1NormExact = 25//32
        L1ErrorValues[i] = abs(testfunc1L1NormExact - WN.norm_l1(testfunc1, xTestGrid, tTestGrid))

        testfunc1L2NormExact = 7//8
        L2ErrorValues[i] = abs(testfunc1L2NormExact - WN.norm_l2(testfunc1, xTestGrid, tTestGrid)^2)
    end
    @test all(d <=0 for d in diff(L1ErrorValues))
    @test all(d <=0 for d in diff(L2ErrorValues))
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

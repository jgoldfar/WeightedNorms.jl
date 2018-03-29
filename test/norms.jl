# norms.jl: tests
function testnormssob()
  include(joinpath(pkgbasedir, "src", "norms.jl"))

  zerofunc(t) = zero(t)
  constfunc(t) = one(t)
  linfunc(t) = t
  sqrfunc(t) = (t * t) / 2
  cubefunc(t) = (t * t * t) / 6
  const testgrid = linspace(0, 1)
  const ntestgrid = length(testgrid)
  const normsjl_epsilon = 1e-10
  @printf "Test grid length %d.\n" ntestgrid
  randgrid = [0.360248873555874,0.9267130090195905,0.6242069066353966,0.6952020993710202,0.6115879542157936,0.3665643627482258,0.5742293978191828,0.21447959790328097,0.012021462291265994,0.28165062937238816,0.10622621222413398,0.3639612507740333,0.29635736283391534,0.6359767569756187,0.1811462873079508,0.0030369911011354223,0.16067902107258947,0.1994895945836861,0.3837955401061268,0.8113715523026992,0.824174763341661,0.1876763556349692,0.18276840303278075,0.14369968274489198,0.880594040860657,0.4814398477649371,0.1296908851814873,0.11372481813930602,0.05334033383759418,0.7609074303660024,0.45130169193572645,0.31846328069047036,0.01764352340109454,0.00862809221975347,0.3436974056023161,0.7986929984281819,0.5694660201293662,0.3812385470909132,0.6263427682936222,0.3921977375206631,0.1332772178126993,0.1459914096177164,0.5041541078140601,0.5327882381461069,0.7586705763719488,0.5892433219015798,0.7339445950652914,0.9655966160777532,0.00033702036775129507,0.3402717836619158,0.07535082269511206,0.50934526901078,0.16191737473788814,0.4190211953285161,0.21201100405864448,0.9909079147894482,0.3803883167081916,0.30943540960962546,0.2577317914775652,0.8004659928355133,0.29541808202411923,0.536470663642928,0.5167528796491578,0.9422365975915121,0.47281680491552525,0.34664708736623995,0.5332253591626657,0.3089401580090052,0.4718255046002755,0.04172552774521221,0.36415352584729344,0.8621131585802342,0.41140600219066803,0.6961709186226412,0.4898439497986047,0.7687979713274884,0.09188723798132092,0.4639372500076777,0.3509960592879109,0.483635026715445,0.5452236688652536,0.5400252275390187,0.7325488065779568,0.9358466736258808,0.5992687849059657,0.9874308789766042,0.2240410825357051,0.14002610429779105,0.007752267386566647,0.07120165485071173,0.008093949177849336,0.8732173611552547,0.5016719336569135,0.06232614806677739,0.4393363101238552,0.13052914834148344,0.05083067651922124,0.35046334846867166,0.08406876904390348,0.06877912313144363]

  function equalitytest()
    for func in sobnormsfunc
      for testfunc in (zerofunc, constfunc, linfunc, sqrfunc, cubefunc, sin, exp)
        if verbose
          @printf "Testing norm %s version equality against data generated from %s.\n" func testfunc
        end
        testdata = map(testfunc, testgrid)
        f1 = func(testdata, testgrid)
        f2 = func(testdata, testgrid[2] - testgrid[1])

        @test isapprox(f1, f2)
        isapprox(f1, f2) || (@printf "\nNonuniform: %f, Uniform: %f, Error: %f\n" f1 f2 (f1 - f2))
      end
      if verbose
        @printf "Testing norm %s version equality against random data.\n" func
      end
      f1 = func(randgrid, testgrid)
      f2 = func(randgrid, testgrid[2] - testgrid[1])
      @test isapprox(f1, f2)
      isapprox(f1, f2) || (@printf "\nNonuniform: %f, Uniform: %f, Error: %f\n" f1 f2 (f1 - f2))
    end #for
    return nothing
  end #function
  @testset "Equality" begin
    equalitytest()
  end

  function knowntest()
    # Norm of zero is zero
    for func in sobnormsfunc
      if verbose
        println("Testing norm ", func, " coercivity.")
      end
      testdata = map(zerofunc, testgrid)
      f1 = func(testdata, testgrid)
      f2 = func(testdata, testgrid[2] - testgrid[1])
      @test isapprox(f1, 0.0)
      @test isapprox(f2, 0.0)
      isapprox(f1, 0.0) || (@printf "\nNonuniform error: %f\n" f1)
      isapprox(f2, 0.0) || (@printf "\nUniform error: %f\n" f2)
      #   Norm is non-negative
      for testfunc in [constfunc, linfunc, sqrfunc, cubefunc, sin, exp]
        if verbose
          @printf "Testing norm %s positivity with function %s.\n" func testfunc
        end
        testdata = map(testfunc, testgrid)
        f1 = func(testdata, testgrid)
        f2 = func(testdata, testgrid[2] - testgrid[1])
        @test f1 >= 0.0
        f1 >= 0.0 || (@printf "\nNonuniform error: %f\n" f1)
        @test f2 >= 0.0 || (@printf "\nUniform error: %f\n" f2)
      end # for testfunc
    end # for sobnormsfunc

    # Other known values
    const known = (
      (constfunc, norm_l2_only, 1.0),
      (constfunc, normsq_w21_only, 0.0),
      (constfunc, normsq_w22_only, 0.0),
      (constfunc, norm_w21_only, 0.0),
      (constfunc, norm_w22_only, 0.0),
      (linfunc, norm_w22_only, 0.0),
      (linfunc, normsq_w22_only, 0.0),
      (linfunc, norm_w21_only, 1.0),
      (sqrfunc, norm_w22_only, 1.0)
      )
    for (testfunc, normfunc, expected) in known
      testdata = map(testfunc, testgrid)
      f1 = normfunc(testdata, testgrid)
      f2 = normfunc(testdata, testgrid[2] - testgrid[1])
      if verbose
        @printf "Testing %s against expected value on %s.\n" normfunc testfunc
      end
      @test (abs(f1 - expected) < normsjl_epsilon)
      (abs(f1 - expected) < normsjl_epsilon) || (@printf "\nNonuniform value: %f, expected: %f, error: %f\n" f1 expected (f1 - expected))
      @test (abs(f2 - expected) < normsjl_epsilon)
      (abs(f2 - expected) < normsjl_epsilon) || (@printf "\nUniform value: %f, expected: %f, error: %f\n" f2 expected (f2 - expected))
    end #for
    return nothing
  end #function
  @testset "Known Values" begin
    knowntest()
  end

  function matrixtest()
    const testmat = ones(ntestgrid, ntestgrid)
    for func in (norm_l2_only, norm_l1, norm_linf)
      if verbose
        @printf "testing %s on matrix input.\n" func
      end
      v1 = func(testmat, testgrid, testgrid)
      @test isapprox(v1, 1.0)
      isapprox(v1, 1.0) || (@printf "\nValue: %f, expected: %f, error: %f\n" v1 1.0 (v1 - 1.0))
    end
    return nothing
  end
  @testset "Matrix Norms" begin
    matrixtest()
  end
  return nothing
end # function
@testset "norms.jl" begin
  testnormssob()
end

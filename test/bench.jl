using WeightedNorms

@static if VERSION < v"0.7-"
    using Compat: range
else
    using Printf
end

## Norm benchmarking
const nsamples = 100
maxNormNameLength = maximum(length(string(s)) for s in WeightedNorms.exportedNormSymbols)
nExportedNorms = length(WeightedNorms.exportedNormSymbols)
for len in [2^n for n in 6:10]
  x = rand(len)
  grid = range(0, stop=1, length=len)
  tau = step(grid)
  tt1 = ones(nExportedNorms)
  tt2 = ones(nExportedNorms)
  for (k, func) in enumerate(WeightedNorms.exportedNorms)
    func(x, grid)
    func(x, tau)
    gc()
    tt1[k] = @elapsed for i in 1:nsamples; func(x, grid); end
    gc()
    tt2[k] = @elapsed for i in 1:nsamples; func(x, tau); end
  end
  for (i, (t1, t2)) in enumerate(zip(tt1, tt2))
      funci = rpad(string(WeightedNorms.exportedNormSymbols[i]), maxNormNameLength, " ")
      @printf "n: %d, func: %s, time 1: %e, time 2: %e, grid is: %f x slower\n" len funci t1 t2 t1/t2
  end
end
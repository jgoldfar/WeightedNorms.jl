const pkgbasedir = dirname(dirname(@__FILE__))
include(joinpath(pkgbasedir, "src", "norms.jl"))

## Norm benchmarking
const nsamples = 10
for len in 2.^(7:10)
  const x = rand(len)
  const grid = linspace(0, 1, len)
  if VERSION >= v"0.4-"
    const tau = step(grid)
  else
    const tau = grid[2] - grid[1]
  end
  for func in sobnormsfunc
    func(x, grid)
    func(x, tau)
    gc()
    timetaken1 = @elapsed for i in 1:nsamples; func(x, grid); end
    gc()
    timetaken2 = @elapsed for i in 1:nsamples; func(x, tau); end
    @printf "n: %d, func: %s, time 1: %e, time 2: %e, ratio: %f\n" len rpad(func, max(map(length, map(string, sobnormsfunc))...), " ") timetaken1 timetaken2 timetaken1/timetaken2
  end
end
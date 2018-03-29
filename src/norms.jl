const honestymode = true # Make true to emulate approximation of Lebesgue norms

##
# Begin general norm routines
##
normsq_l2{T<:Real}(x::Vector{T}) = sumabs2(x)
normsq_l2{T<:Real}(x::Matrix{T}) = sumabs2(collect(x))
diffsq_l2{T<:Real}(x1::Matrix{T}, x2::Matrix{T}) = sumabs2(collect(x1 - x2))
norm_l2{T<:Real}(x::VecOrMat{T}) = sqrt(normsq_l2(x))
diff_l2{T<:Real}(x1::Matrix{T}, x2::Matrix{T}) = sqrt(diffsq_l2(x1, x2))
norm_linf{T<:Real}(x::VecOrMat{T}) = vecnorm(x, Inf)
diff_linf{T<:Real}(x1::Vector{T}, x2::Vector{T}) = diff_linf(x1 - x2)
diff_linf{T<:Real}(x1::Matrix{T}, x2::Matrix{T}) = diff_linf(x1 - x2)
norm_linf{T<:Real}(x::Vector{T}, grid::Vector) = maxabs(x)
norm_linf{T<:Real}(x::Matrix{T}, grid1, grid2) = vecnorm(x, Inf)
diff_linf{T<:Real}(x1::Matrix{T}, x2::Matrix{T},
                   grid1, grid2) = norm_linf(x1 - x2, grid1, grid2)

function normsq_l2{T<:Real}(x::Matrix{T}, grid1::Vector{T}, grid2::Vector{T})
  xnext, v0 = grid2[1], zero(T)
  const sz1, sz2 = size(x)
  for j in 2:sz2
    xcurr, xnext = xnext, grid2[j]
    dx = xnext - xcurr
    tnext = grid1[1]
    for i in 2:sz1
      tcurr, tnext = tnext, grid1[i]
      v0 += abs2(x[i, j]) * dx * (tnext - tcurr)
    end
  end
  v0
end
function normsq_l2{T<:Real}(x::Matrix{T}, dx::Real, dt::Real)
  v0 = zero(T)
  const sz1, sz2 = size(x)
  for j in 2:sz2
    for i in 2:sz1
      v0 += abs2(x[i, j])
    end
  end
  v0 * dt * dx
end
if VERSION >= v"0.4-"
  normsq_l2{T<:Real}(x::Matrix{T}, grid1::LinSpace, grid2::LinSpace) = normsq_l2(x, step(grid1), step(grid2)) 
end
norm_l2_only{T<:Real}(x::Matrix{T}, grid1, grid2) = sqrt(normsq_l2(x, grid1, grid2))

function diffsq_l2{T<:Real}(x1::Matrix{T}, x2::Matrix{T}, grid1::Vector{T}, grid2::Vector{T})
  xnext, v0 = grid2[1], zero(T)
  const sz1, sz2 = size(x1)
  for j in 2:sz2
    xcurr, xnext = xnext, grid2[j]
    dx = xnext - xcurr
    tnext = grid1[1]
    for i in 2:sz1
      tcurr, tnext = tnext, grid1[i]
      @inbounds v0 += abs2(x1[i, j] - x2[i, j]) * dx * (tnext - tcurr)
    end
  end
  v0
end
function diffsq_l2{T<:Real}(x1::Matrix{T}, x2::Matrix{T}, dx::Real, dt::Real)
  v0 = zero(T)
  const sz1, sz2 = size(x1)
  for j in 2:sz2
    for i in 2:sz1
      @inbounds v0 += abs2(x1[i, j] - x2[i, j])
    end
  end
  v0 * dx * dt
end
if VERSION >= v"0.4-"
  diffsq_l2{T<:Real}(x1::Matrix{T}, x2::Matrix{T}, grid1::LinSpace, grid2::LinSpace) = diffsq_l2(x1, x2, step(grid1), step(grid2))
end

diff_l2_only{T<:Real}(x1::Matrix{T}, x2::Matrix{T}, grid1::Vector, grid2::Vector) = sqrt(diffsq_l2(x1, x2, grid1, grid2))

function norm_l1{T<:Real}(x::Matrix{T}, grid1::Vector{T}, grid2::Vector{T})
  xnext, v0 = grid2[1], zero(T)
  const sz1, sz2 = size(x)
  for j in 2:sz2
    xcurr, xnext = xnext, grid2[j]
    dx = xnext - xcurr
    tnext = grid1[1]
    for i in 2:sz1
      tcurr, tnext = tnext, grid1[i]
      @inbounds v0 += abs(x[i, j]) * dx * (tnext - tcurr)
    end
  end
  return v0
end
function norm_l1{T<:Real}(x::Matrix{T}, dx::Real, dt::Real)
  v0 = zero(T)
  const sz1, sz2 = size(x)
  for j in 2:sz2
    for i in 2:sz1
      @inbounds v0 += abs(x[i, j])
    end
  end
  return v0 * dx * dt
end
if VERSION >= v"0.4-"
  norm_l1{T<:Real}(x::Matrix{T}, grid1::LinSpace, grid2::LinSpace) = norm_l1(x, step(grid1), step(grid2))
end

function diff_l1{T<:Real}(x1::Matrix{T}, x2::Matrix{T}, grid1::Vector{T}, grid2::Vector{T})
  xnext, v0 = grid2[1], zero(T)
  const sz1, sz2 = size(x1)
  for j in 2:sz2
    xcurr, xnext = xnext, grid2[j]
    dx = xnext - xcurr
    tnext = grid1[1]
    for i in 2:sz1
      tcurr, tnext = tnext, grid1[i]
      @inbounds v0 += abs(x1[i, j] - x2[i, j]) * dx * (tnext - tcurr)
    end
  end
  v0
end
function diff_l1{T<:Real}(x1::Matrix{T}, x2::Matrix{T}, dx::Real, dt::Real)
  v0 = zero(T)
  const sz1, sz2 = size(x1)
  for j in 2:sz2
    for i in 2:sz1
      @inbounds v0 += abs(x1[i, j] - x2[i, j])
    end
  end
  v0 * dx * dt
end
if VERSION >= v"0.4-"
  diff_l1{T<:Real}(x1::Matrix{T}, x2::Matrix{T}, grid1::LinSpace, grid2::LinSpace) = diff_l1(x1, x2, step(grid1), step(grid2))
end
##
# End general norm routines
##


export norm_l2_only, norm_w21_only, norm_w22_only, norm_w21, norm_w22, norm_null
if honestymode
  function normsq_l2{T<:Real}(x::Vector{T}, grid::Vector{T})
    gnext, v0 = grid[1], zero(T)
    for i in 2:length(x)
      gcurr, gnext = gnext, grid[i]
      v0 += abs2(x[i]) * (gnext - gcurr)
    end
    v0
  end
  function diffsq_l2{T<:Real}(x1::Vector{T}, x2::Vector{T}, grid::Vector{T})
    gnext, v0 = grid[1], zero(T)
    for i in 2:length(x1)
      gcurr, gnext = gnext, grid[i]
      @inbounds v0 += abs2(x1[i] - x2[i]) * (gnext - gcurr)
    end
    v0
  end
  function normsq_l2{T<:Real}(x::Vector{T}, dx::Real)
    v0 = zero(T)
    for i in 2:length(x)
      @inbounds v0 += abs2(x[i])
    end
    v0 * dx
  end
  function diffsq_l2{T<:Real}(x1::Vector{T}, x2::Vector{T}, dx::Real)
    v0 = zero(T)
    for i in 2:length(x1)
      @inbounds v0 += abs2(x1[i] - x2[i])
    end
    v0 * dx
  end
else
  function normsq_l2{T<:Real}(x::Vector{T}, grid::Vector{T})
    gnext = grid[1]
    @inbounds v0 = abs2(x[1]) * (grid[2] - gnext)
    for i in 2:length(x)
      gcurr, gnext = gnext, grid[i]
      v0 += abs2(x[i]) * (gnext - gcurr)
    end
    v0
  end
  function diffsq_l2{T<:Real}(x1::Vector{T}, x2::Vector{T}, grid::Vector{T})
    gnext = grid[1]
    @inbounds v0 = abs2(x2[1] - x1[1]) * (grid[2] - gnext)
    for i in 2:length(x1)
      gcurr, gnext = gnext, grid[i]
      @inbounds v0 += abs2(x1[i] - x2[i]) * (gnext - gcurr)
    end
    v0
  end
  function normsq_l2{T<:Real}(x::Vector{T}, dx::Real)
    v0 = abs2(x[1])
    for i in 2:length(x)
      @inbounds v0 += abs2(x[i])
    end
    v0 * dx
  end
  function diffsq_l2{T<:Real}(x1::Vector{T}, x2::Vector{T}, dx::Real)
    v0 = abs2(x1[1] - x2[1])
    for i in 2:length(x1)
      @inbounds v0 += abs2(x1[i] - x2[i])
    end
    v0 * dx
  end
end
norm_l2_only{T<:Real}(x::Vector{T}, grid::Vector{T}) = sqrt(normsq_l2(x, grid))
norm_l2_only{T<:Real}(x::Vector{T}, dx::Real) = sqrt(normsq_l2(x, dx))
if VERSION >= v"0.4-"
  norm_l2_only{T<:Real}(x::Vector{T}, grid::LinSpace) = norm_l2_only(x, step(grid))
end

diff_l2{T<:Real}(x1::Vector{T}, x2::Vector{T}, grid::Vector{T}) = sqrt(diffsq_l2(x1, x2, grid))
diff_l2{T<:Real}(x1::Vector{T}, x2::Vector{T}, dx::Real) = sqrt(diffsq_l2(x1, x2, dt))
if VERSION >= v"0.4-"
  diff_l2{T<:Real}(x1::Vector{T}, x2::Vector{T}, grid::LinSpace) = sqrt(diffsq_l2(x1, x2, step(grid)))
end

##
# Begin discrete Sobolev norm routines
##

function normsq_w21_only{T<:Real}(x::Vector{T}, grid::Vector{T})
  xnext, gnext, v0 = x[1], grid[1], zero(T)
  for i in 2:length(x)
    xcurr, xnext, gcurr, gnext = xnext, x[i], gnext, grid[i]
    v0 += abs2(xnext - xcurr) / (gnext - gcurr)
  end
  v0
end
function normsq_w21_only{T<:Real}(x::Vector{T}, dx::Real) # Much faster for large (i.e. >512 points.)
  xnext, v0 = x[1], zero(T)
  for i in 2:length(x)
    xcurr, xnext = xnext, x[i]
    v0 += abs2(xnext - xcurr)
  end
  v0 / dx
end
if VERSION >= v"0.4-"
  normsq_w21_only{T<:Real}(x::Vector{T}, grid::LinSpace) = normsq_w21_only(x, step(grid))
end

norm_w21_only{T<:Real}(x::Vector{T}, grid::Vector{T}) = sqrt(normsq_w21_only(x, grid))
norm_w21_only{T<:Real}(x::Vector{T}, dx::Real) = sqrt(normsq_w21_only(x, dx))
if VERSION >= v"0.4-"
  norm_w21_only{T<:Real}(x::Vector{T}, grid::LinSpace) = sqrt(normsq_w21_only(x, step(grid)))
end

normsq_w21{T<:Real}(x::Vector{T}, grid::Vector{T}) = normsq_l2(x, grid) + normsq_w21_only(x, grid)
normsq_w21{T<:Real}(x::Vector{T}, dx::Real) = normsq_l2(x, dx) + normsq_w21_only(x, dx)
if VERSION >= v"0.4-"
  normsq_w21{T<:Real}(x::Vector{T}, grid::LinSpace) = normsq_l2(x, step(grid)) + normsq_w21_only(x, step(grid))
end

norm_w21{T<:Real}(x::Vector{T}, grid::Vector{T}) = sqrt(normsq_w21(x, grid))
norm_w21{T<:Real}(x::Vector{T}, dx::Real) = sqrt(normsq_w21(x, dx))
if VERSION >= v"0.4-"
  norm_w21{T<:Real}(x::Vector{T}, grid::LinSpace) = sqrt(normsq_w21(x, step(grid)))
end

function normsq_w22_only{T<:Real}(x::Vector{T}, grid::Vector{T})
  #   Note: to handle the general case with full accuracy, a more involved algorithm
  #   is required, leading this to be quite a bit more computationally intensive...
  const nx = length(x)
#   v0 = zero(T)
  if nx < 3
    return NaN
#     return v0
  end
  xprev, xcurr, xnext, dx1, dx2 = x[1], x[2], x[3], grid[2] - grid[1], grid[3] - grid[2]
  rdx1, rdx2, rdxs = one(T) / dx1, one(T) / dx2, one(T)/(dx1 + dx2)
  v0 = 4abs2(rdx1 * rdxs * xprev - rdx1 * rdx2 * xcurr + rdx2 * rdxs * xnext)
  if nx == 3
    return v0 * dx1
  end
  for i in 3:(nx - 1)
    xprev, xcurr, xnext, dx1, dx2 = xcurr, xnext, x[i + 1], dx2, grid[i + 1] - grid[i]
    rdx1, rdx2, rdxs = rdx2, one(T) / dx2, one(T)/(dx1 + dx2)
    v0 += 4*abs2(rdx1 * rdxs * xprev - rdx1 * rdx2 * xcurr + rdx2 * rdxs * xnext)
  end
  v0 * dx1
end
function normsq_w22_only{T<:Real}(x::Vector{T}, dx::Real)
  const nx = length(x)
#   v0 = zero(T)
  if nx < 3
    return NaN
#     return v0
  end
  xprev, xcurr, xnext = x[1], x[2], x[3]
  v0 = abs2(xprev - 2xcurr + xnext)
  if nx == 3
    return v0 / (dx * dx * dx)
  end
  for i in 3:(nx - 1)
    xprev, xcurr, xnext = xcurr, xnext, x[i + 1]
    v0 += abs2(xprev - 2xcurr + xnext)
  end
  v0 / (dx * dx * dx)
end
if VERSION >= v"0.4-"
  normsq_w22_only{T<:Real}(x::Vector{T}, grid::LinSpace) = normsq_w22_only(x, step(grid))
end

norm_w22_only{T<:Real}(x::Vector{T}, grid::Vector{T}) = sqrt(normsq_w22_only(x, grid))
norm_w22_only{T<:Real}(x::Vector{T}, dx::Real) = sqrt(normsq_w22_only(x, dx))

if VERSION >= v"0.4-"
  norm_w22_only{T<:Real}(x::Vector{T}, grid::LinSpace) = sqrt(normsq_w22_only(x, step(grid)))
end

normsq_w22{T<:Real}(x::Vector{T}, grid::Vector{T}) = normsq_w21(x, grid) + normsq_w22_only(x, grid)
normsq_w22{T<:Real}(x::Vector{T}, dx::Real) = normsq_w21(x, dx) + normsq_w22_only(x, dx)

if VERSION >= v"0.4-"
  normsq_w22{T<:Real}(x::Vector{T}, grid::LinSpace) = normsq_w21(x, step(grid)) + normsq_w22_only(x, step(grid))
end

norm_w22{T<:Real}(x::Vector{T}, dx::Real) = sqrt(normsq_w22(x, dx))
norm_w22{T<:Real}(x::Vector{T}, grid::Vector{T}) = sqrt(normsq_w22(x, grid))
if VERSION >= v"0.4-"
  norm_w22{T<:Real}(x::Vector{T}, grid::LinSpace) = sqrt(normsq_w22(x, step(grid)))
end

norm_null{T<:Real}(x::Vector{T}, grid::Vector{T}) = zero(T)
norm_null{T<:Real}(x::Vector{T}, dx::Real) = zero(T)
if VERSION >= v"0.4-"
  norm_null{T<:Real}(x::Vector{T}, grid::LinSpace) = zero(T)
end

# const sobnorms = [:norm_l2_only, :norm_w21, :norm_w22, :norm_w21_only, :norm_w22_only, :norm_null]
const sobnormsfunc = (norm_l2_only, norm_w21, norm_w22, norm_w21_only, norm_w22_only, norm_null)

##
# End discrete Sobolev norm routines
##
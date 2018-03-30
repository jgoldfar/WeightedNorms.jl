##
# Begin general norm routines
##
function normsq_l2{T<:Real}(x::Matrix{T}, grid1::Vector{T}, grid2::Vector{T})
    xnext, v0 = grid2[1], zero(T)
    sz1, sz2 = size(x)
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
    sz1, sz2 = size(x)
    for j in 2:sz2
        for i in 2:sz1
            v0 += abs2(x[i, j])
        end
    end
    v0 * dt * dx
end
normsq_l2(x::Matrix, grid1::AbstractRange, grid2::AbstractRange) = normsq_l2(x, step(grid1), step(grid2))
norm_l2_only(x::Matrix, grid1, grid2) = sqrt(normsq_l2(x, grid1, grid2))

function norm_l1{T<:Real}(x::Matrix{T}, grid1::Vector{T}, grid2::Vector{T})
    xnext, v0 = grid2[1], zero(T)
    sz1, sz2 = size(x)
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
    sz1, sz2 = size(x)
    for j in 2:sz2
        for i in 2:sz1
            @inbounds v0 += abs(x[i, j])
        end
    end
    return v0 * dx * dt
end
norm_l1{T<:Real}(x::Matrix{T}, grid1::AbstractRange, grid2::AbstractRange) = norm_l1(x, step(grid1), step(grid2))
##
# End general norm routines
##


export norm_l2_only, norm_w21_only, norm_w22_only, norm_w21, norm_w22, norm_null

function normsq_l2{T<:Real}(x::Vector{T}, grid::Vector{T})
    gnext = grid[1]
    v0 = zero(T)
    for i in 2:length(x)
        gcurr, gnext = gnext, grid[i]
        v0 += abs2(x[i]) * (gnext - gcurr)
    end
    v0
end
function normsq_l2{T<:Real}(x::Vector{T}, dx::Real)
    v0 = zero(T)
    for i in 2:length(x)
        v0 += abs2(x[i])
    end
    v0 * dx
end
norm_l2_only(x::Vector, grid::Vector) = sqrt(normsq_l2(x, grid))
norm_l2_only(x::Vector, dx::Real) = sqrt(normsq_l2(x, dx))

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

norm_w21_only(x::Vector, grid::Vector) = sqrt(normsq_w21_only(x, grid))
norm_w21_only(x::Vector, dx::Real) = sqrt(normsq_w21_only(x, dx))

normsq_w21(x::Vector, grid::Vector) = normsq_l2(x, grid) + normsq_w21_only(x, grid)
normsq_w21(x::Vector, dx::Real) = normsq_l2(x, dx) + normsq_w21_only(x, dx)

norm_w21(x::Vector, grid::Vector) = sqrt(normsq_w21(x, grid))
norm_w21(x::Vector, dx::Real) = sqrt(normsq_w21(x, dx))

function normsq_w22_only{T<:Real}(x::Vector{T}, grid::Vector{T})
    #   Note: to handle the general case with full accuracy, a more involved algorithm
    #   is required, leading this to be quite a bit more computationally intensive...
    nx = length(x)
    #   v0 = zero(T)
    if nx < 3
        return NaN
        #     return v0
    end
    xprev, xcurr, xnext, dx1, dx2 = x[1], x[2], x[3], grid[2] - grid[1], grid[3] - grid[2]
    rdx1, rdx2, rdxs = 1 / dx1, 1 / dx2, 1/(dx1 + dx2)
    v0 = 4*abs2(rdx1 * rdxs * xprev - rdx1 * rdx2 * xcurr + rdx2 * rdxs * xnext)
    if nx == 3
        return v0 * dx1
    end
    for i in 3:(nx - 1)
        xprev, xcurr, xnext, dx1, dx2 = xcurr, xnext, x[i + 1], dx2, grid[i + 1] - grid[i]
        rdx1, rdx2, rdxs = rdx2, 1 / dx2, 1 / (dx1 + dx2)
        v0 += 4*abs2(rdx1 * rdxs * xprev - rdx1 * rdx2 * xcurr + rdx2 * rdxs * xnext)
    end
    v0 * dx1
end
function normsq_w22_only(x::Vector, dx::Real)
    nx = length(x)
    if nx < 3
        return NaN
    end
    xprev, xcurr, xnext = x[1], x[2], x[3]
    v0 = abs2(xprev - 2*xcurr + xnext)
    if nx == 3
        return v0 / (dx * dx * dx)
    end
    for i in 3:(nx - 1)
        xprev, xcurr, xnext = xcurr, xnext, x[i + 1]
        v0 += abs2(xprev - 2*xcurr + xnext)
    end
    v0 / (dx * dx * dx)
end

norm_w22_only(x::Vector, grid::Vector) = sqrt(normsq_w22_only(x, grid))
norm_w22_only(x::Vector, dx::Real) = sqrt(normsq_w22_only(x, dx))

normsq_w22(x::Vector, grid::Vector) = normsq_w21(x, grid) + normsq_w22_only(x, grid)
normsq_w22(x::Vector, dx::Real) = normsq_w21(x, dx) + normsq_w22_only(x, dx)

norm_w22(x::Vector, grid::Vector) = sqrt(normsq_w22(x, grid))
norm_w22(x::Vector, dx::Real) = sqrt(normsq_w22(x, dx))

norm_null{T<:Real}(x::Vector{T}, ::ANY) = zero(T)

const exportedNorms = [norm_l2_only,
                       normsq_w21, norm_w21, normsq_w21_only, norm_w21_only,
                       normsq_w22, norm_w22, normsq_w22_only, norm_w22_only,
                       norm_null]
const exportedNormSymbols = [:norm_l2_only,
                             :normsq_w21, :norm_w21, :normsq_w21_only, :norm_w21_only,
                             :normsq_w22, :norm_w22, :normsq_w22_only, :norm_w22_only,
                             #:norm_null
                             ]
##
# End discrete Sobolev norm routines
##

## Define functions on AbstractRanges using definitions with constant grid spacing
for fn in exportedNormSymbols
    @eval ($fn)(x::Vector, grid::AbstractRange) = ($fn)(x, step(grid))
end

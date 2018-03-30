## Norms
function space_trace_norm{T<:Real}(u1::Matrix{T}, xgrid::Vector, Nspacep::Int, ind::Int)
    out = zero(T)
    for i in 2:Nspacep
        out += (xgrid[i] - xgrid[i - 1]) * abs2(u1[i, ind])
    end
    sqrt(out)
end
function space_trace_norm{T<:Real}(u1::Matrix{T}, xgrid::AbstractRange, Nspacep::Int, ind::Int)
    out = zero(T)
    for i in 2:Nspacep
        out += abs2(u1[i, ind])
    end
    sqrt(step(xgrid) * out)
end

function space_trace_dist{T<:Real}(u1::Matrix{T}, v::Vector{T}, xgrid::Vector, Nspacep::Int, ind::Int)
    out = zero(T)
    for i in 2:Nspacep
        out += (xgrid[i] - xgrid[i - 1]) * abs2(u1[i, ind] - v[i])
    end
    sqrt(out)
end
function space_trace_dist{T<:Real}(u1::Matrix{T}, v::Vector{T}, xgrid::AbstractRange, Nspacep::Int, ind::Int)
    out = zero(T)
    for i in 2:Nspacep
        out += abs2(u1[i, ind] - v[i])
    end
    sqrt(step(xgrid) * out)
end

function space_dist{T<:Real}(u1::Matrix{T}, u2::Matrix{T}, xgrid::Vector, k::Real)
    NX, NT = size(u1)
    dout = zero(T)
    for i in 1:NT
        for j in 2:NX
            dout += (xgrid[j] - xgrid[j - 1]) * abs2(u1[j, i] - u2[j, i])
        end
    end
    sqrt(k * dout)
end
function space_dist{T<:Real}(u1::Matrix{T}, u2::Matrix{T}, xgrid::AbstractRange, k::Real)
    NX, NT = size(u1)
    dout = zero(T)
    for i in 1:NT
        for j in 2:NX
            dout += abs2(u1[j, i] - u2[j, i])
        end
    end
    sqrt(k * step(xgrid) * dout)
end

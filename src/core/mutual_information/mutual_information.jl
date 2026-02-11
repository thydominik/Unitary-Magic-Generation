"""
    mutual_information

Mutual information utilities.

This module provides simple, dependency-light mutual information estimators that are useful
for quick analysis in this repository.

All identifiers and docstrings are ASCII-only.

Implemented
- Discrete mutual information given a joint probability table.
- Histogram-based mutual information estimator for 1D samples (plug-in estimator).

Notes
- All results are in bits (log base 2).
- Histogram estimation is sensitive to binning; treat results as qualitative unless you validate.
"""
module mutual_information

using Random

export mutual_information_discrete,
    mutual_information_histogram

"""
    mutual_information_discrete(pxy) -> Float64

Compute mutual information I(X;Y) in bits from a joint probability table.

Parameters
- `pxy`: A nonnegative matrix with sum(pxy) approximately 1.

Returns
- Mutual information in bits.

Definition
I(X;Y) = sum_{i,j} pxy[i,j] * log2( pxy[i,j] / (px[i] * py[j]) )
where px = sum over columns, py = sum over rows.
"""
function mutual_information_discrete(pxy::AbstractMatrix{<:Real})::Float64
    any(pxy .< 0) && throw(ArgumentError("pxy must be nonnegative"))

    s = float(sum(pxy))
    s > 0 || return 0.0

    p = float.(pxy) ./ s
    px = vec(sum(p, dims=2))
    py = vec(sum(p, dims=1))

    mi = 0.0
    @inbounds for i in axes(p, 1)
        for j in axes(p, 2)
            pij = p[i, j]
            if pij > 0
                denom = px[i] * py[j]
                if denom > 0
                    mi += pij * log2(pij / denom)
                end
            end
        end
    end

    return mi
end

"""Return bin edges for `n_bins` bins spanning [a,b]."""
function _edges(a::Real, b::Real, n_bins::Integer)::Vector{Float64}
    n_bins >= 1 || throw(ArgumentError("n_bins must be >= 1"))
    a_f = Float64(a)
    b_f = Float64(b)
    a_f == b_f && (b_f = a_f + 1.0)
    return collect(range(a_f, b_f; length=n_bins + 1))
end

"""Map value to bin index in 1..n_bins, or return 0 if outside [edges[1], edges[end])."""
function _bin_index(val::Real, edges::AbstractVector{<:Real})::Int
    v = Float64(val)
    if v < edges[1] || v > edges[end]
        return 0
    end
    # Rightmost edge is inclusive; shift it to the last bin.
    if v == edges[end]
        return length(edges) - 1
    end
    k = searchsortedlast(edges, v)
    return max(1, min(k, length(edges) - 1))
end

"""
    mutual_information_histogram(x, y; n_bins=32, x_min=nothing, x_max=nothing, y_min=nothing, y_max=nothing) -> Float64

Estimate mutual information I(X;Y) in bits from samples using a 2D histogram plug-in estimator.

Parameters
- `x`, `y`: Equal-length vectors of real samples.
- `n_bins`: Number of bins per axis.
- `x_min`, `x_max`, `y_min`, `y_max`: Optional fixed ranges. If not provided, min/max from data are used.

Returns
- Mutual information in bits.

Notes
- Samples falling outside the specified ranges are ignored.
- This estimator is biased for finite samples.
"""
function mutual_information_histogram(
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real};
    n_bins::Integer=32,
    x_min=nothing,
    x_max=nothing,
    y_min=nothing,
    y_max=nothing,
)::Float64
    length(x) == length(y) || throw(DimensionMismatch("x and y must have the same length"))
    n = length(x)
    n > 0 || return 0.0

    xmin = x_min === nothing ? minimum(x) : x_min
    xmax = x_max === nothing ? maximum(x) : x_max
    ymin = y_min === nothing ? minimum(y) : y_min
    ymax = y_max === nothing ? maximum(y) : y_max

    ex = _edges(xmin, xmax, n_bins)
    ey = _edges(ymin, ymax, n_bins)

    counts = zeros(Float64, n_bins, n_bins)
    used = 0

    @inbounds for k in 1:n
        ix = _bin_index(x[k], ex)
        iy = _bin_index(y[k], ey)
        if ix != 0 && iy != 0
            counts[ix, iy] += 1.0
            used += 1
        end
    end

    used > 0 || return 0.0

    pxy = counts ./ used
    return mutual_information_discrete(pxy)
end

end # module mutual_information

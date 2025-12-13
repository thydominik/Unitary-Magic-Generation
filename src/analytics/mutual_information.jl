"""
    MutualInformationAnalysis

Module for computing mutual information between two distributions using histogram-based estimation.
Provides both fixed-range and normalized (data-range) variants.
"""
module MutualInformationAnalysis

using StatsBase: fit, Histogram
using ..NumericalIntegration: double_integral_trapz

export MI, MIn

"""
    MI(x::Vector, y::Vector, N_bins::Int)::Tuple

Compute mutual information using fixed range [0, 10] with trapezoidal integration.

# Arguments
- `x::Vector`: First data array
- `y::Vector`: Second data array  
- `N_bins::Int`: Number of bins for histogram

# Returns
- Tuple containing:
  - `mi::Float64`: Mutual information value
  - `px_individual::Vector`: Marginal histogram of x
  - `px_marginal::Vector`: Marginal probabilities of x
  - `py_individual::Vector`: Marginal histogram of y
  - `py_marginal::Vector`: Marginal probabilities of y
  - `pxy::Array{Float64,2}`: Joint probability distribution
  - `Integrand::Array{Float64,2}`: The integrand used in MI computation

# Theory
MI is computed as:
```
MI(X;Y) = ∫∫ p(x,y) log₂(p(x,y) / (p(x)p(y))) dx dy
```

# Notes
- Uses fixed range [0, 10] for both variables
- Probability density is normalized by bin widths
- Zero-valued bins are skipped to avoid log(0)
"""
function MI(x::Vector, y::Vector, N_bins::Int)
    m = length(x)

    # Define fixed range bins
    bins_x = range(0, 10, N_bins)
    d_bins_x = bins_x[2] - bins_x[1]

    bins_y = range(0, 10, N_bins)
    d_bins_y = bins_y[2] - bins_y[1]

    # Compute joint histogram and normalize by density
    joint_hist = fit(Histogram, (x, y), (bins_x, bins_y))
    pxy = joint_hist.weights ./ (m * d_bins_x * d_bins_y)

    # Marginal probabilities for x
    px_marginal = sum(pxy, dims=2)
    h_x = fit(Histogram, x, bins_x)
    px_individual = h_x.weights / (m * d_bins_x)

    # Marginal probabilities for y
    py_marginal = sum(pxy, dims=1)
    h_y = fit(Histogram, y, bins_y)
    py_individual = h_y.weights / (m * d_bins_y)

    # Compute integrand: p(x,y) * log₂(p(x,y) / (p(x)p(y)))
    Integrand = zeros(length(px_individual), length(py_individual))
    for i in 1:length(px_individual)
        for j in 1:length(py_individual)
            if pxy[i, j] > 0
                Integrand[i, j] = pxy[i, j] * log2(pxy[i, j] / (px_individual[i] * py_individual[j]))
            end
        end
    end
    
    mi = double_integral_trapz(Integrand, bins_x[1:end-1], bins_y[1:end-1])

    return mi, px_individual, px_marginal, py_individual, py_marginal, pxy, Integrand
end


"""
    MIn(x::Vector, y::Vector, N_bins::Int)::Tuple

Compute mutual information using data-driven range (normalized MI).
Bins are determined by actual data range rather than fixed range.

# Arguments
- `x::Vector`: First data array
- `y::Vector`: Second data array
- `N_bins::Int`: Number of bins for histogram

# Returns
- Tuple containing same as MI():
  - `mi::Float64`: Mutual information value
  - `px_individual::Vector`: Marginal histogram of x
  - `px_marginal::Vector`: Marginal probabilities of x  
  - `py_individual::Vector`: Marginal histogram of y
  - `py_marginal::Vector`: Marginal probabilities of y
  - `pxy::Array{Float64,2}`: Joint probability distribution
  - `Integrand::Array{Float64,2}`: The integrand used in MI computation

# Notes
- Uses data range [min(x), max(x)] and [min(y), max(y)]
- More suitable for data with non-uniform support
- Normalization includes factor of m for discrete sample correction
"""
function MIn(x::Vector, y::Vector, N_bins::Int)
    m = length(x)

    # Data-driven bin ranges
    bins_x = range(minimum(x), maximum(x), N_bins)
    d_bins_x = bins_x[2] - bins_x[1]

    bins_y = range(minimum(y), maximum(y), N_bins)
    d_bins_y = bins_y[2] - bins_y[1]

    # Compute joint histogram and normalize by density
    joint_hist = fit(Histogram, (x, y), (bins_x, bins_y))
    pxy = joint_hist.weights ./ (m * d_bins_x * d_bins_y)

    # Marginal probabilities for x
    px_marginal = sum(pxy, dims=2)
    h_x = fit(Histogram, x, bins_x)
    px_individual = h_x.weights / (m * d_bins_x)

    # Marginal probabilities for y
    py_marginal = sum(pxy, dims=1)
    h_y = fit(Histogram, y, bins_y)
    py_individual = h_y.weights / (m * d_bins_y)

    # Compute integrand with sample count normalization
    Integrand = zeros(length(px_individual), length(py_individual))
    for i in 1:length(px_individual)
        for j in 1:length(py_individual)
            if pxy[i, j] > 0
                Integrand[i, j] = pxy[i, j] / m * log2(pxy[i, j] * m / (px_individual[i] * py_individual[j]))
            end
        end
    end
    
    mi = double_integral_trapz(Integrand, bins_x[1:end-1], bins_y[1:end-1])

    return mi, px_individual, px_marginal, py_individual, py_marginal, pxy, Integrand
end

end  # module MutualInformationAnalysis

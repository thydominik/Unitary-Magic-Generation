"""
    MutualInformationAnalysis

Module for computing mutual information between two distributions using histogram-based 
estimation.

This module provides statistical analysis tools for quantifying dependence between variables
through information-theoretic measures. Two variants are provided:

1. **MI**: Fixed-range mutual information (assumes data in [0, 10])
2. **MIn**: Data-driven mutual information (normalizes to actual data range)

# Key Features

- Histogram-based probability estimation from discrete samples
- 2D numerical integration using trapezoidal rule (from NumericalIntegration module)
- Proper normalization by bin widths for probability density
- Handles singular/zero probability cases to avoid log(0) errors
- Returns complete probability distributions for analysis

# Exported Functions

- `MI`: Compute mutual information with fixed range [0, 10]
- `MIn`: Compute mutual information with data-driven range

# Theory Background

Mutual Information (MI) quantifies the amount of information that one random variable 
carries about another.

## Mathematical Definition

```
MI(X;Y) = ∫∫ p(x,y) log₂(p(x,y) / (p(x)p(y))) dx dy
       = H(X) + H(Y) - H(X,Y)
```

where:
- p(x,y) is the joint probability distribution
- p(x) and p(y) are marginal distributions
- H(X), H(Y), H(X,Y) are Shannon entropies

## Interpretation

- MI = 0: X and Y are independent
- MI > 0: X and Y are dependent
- MI increases with stronger dependence
- Measured in bits when using log₂

# Implementation Details

## Histogram-Based Estimation

Since we only have samples, we estimate p(x,y) using histograms:

1. Create 2D histogram of (x,y) samples into bins
2. Normalize histogram counts to get probability densities
3. Compute marginals by summing over dimensions
4. Compute integrand = p(x,y) * log₂(p(x,y) / (p(x)p(y)))
5. Integrate using 2D trapezoidal rule

## Probability Density Normalization

```
p_density(i,j) = histogram_count(i,j) / (N_samples * bin_width_x * bin_width_y)
```

This ensures the integral of p(x,y) over the domain equals 1.
"""
module MutualInformationAnalysis

# ============================================================================
# MODULE IMPORTS
# ============================================================================
# External dependencies
using StatsBase: fit, Histogram
# Internal dependencies - NOTE: calling functions FROM this module!
using ..NumericalIntegration: double_integral_trapz

# ============================================================================
# EXPORTS
# ============================================================================
# Public API functions that users should call
export MI, MIn

# ============================================================================
# PRIVATE HELPER FUNCTIONS
# ============================================================================
# These are internal functions used by MI and MIn

"""
    _compute_histogram_bins(
        data::AbstractVector, 
        n_bins::Int, 
        fixed_range::Union{Tuple{Real,Real}, Nothing}=nothing
    )::Tuple{AbstractRange, Float64}

Internal helper: compute histogram bins and bin width.

# Arguments
- `data::AbstractVector`: Data to bin
- `n_bins::Int`: Number of bins to create
- `fixed_range::Union{Tuple{Real,Real}, Nothing}`: Optional fixed range (min, max)

# Returns
- Tuple of (bins, bin_width)
"""
function _compute_histogram_bins(
    data::AbstractVector, 
    n_bins::Int, 
    fixed_range::Union{Tuple{Real,Real}, Nothing}=nothing
)::Tuple{AbstractRange, Float64}
    
    if fixed_range !== nothing
        # Use provided fixed range
        bins = range(fixed_range[1], fixed_range[2], n_bins)
    else
        # Use data-driven range
        bins = range(minimum(data), maximum(data), n_bins)
    end
    
    # Compute bin width (assumed uniform)
    bin_width = bins[2] - bins[1]
    
    return bins, bin_width
end

"""
    _compute_probability_distributions(
        x::AbstractVector,
        y::AbstractVector,
        bins_x::AbstractRange,
        bins_y::AbstractRange
    )::Tuple{Matrix, Vector, Vector, Vector, Vector}

Internal helper: compute joint and marginal probability distributions.

# Returns
- Tuple of (pxy, px_individual, px_marginal, py_individual, py_marginal)
"""
function _compute_probability_distributions(
    x::AbstractVector,
    y::AbstractVector,
    bins_x::AbstractRange,
    bins_y::AbstractRange
)::Tuple{Matrix, Vector, Vector, Vector, Vector}
    
    m = length(x)
    d_bins_x = bins_x[2] - bins_x[1]
    d_bins_y = bins_y[2] - bins_y[1]
    
    # Create 2D histogram and normalize to probability density
    joint_hist = fit(Histogram, (x, y), (bins_x, bins_y))
    pxy = joint_hist.weights ./ (m * d_bins_x * d_bins_y)
    
    # Marginal probabilities for X
    px_marginal = sum(pxy, dims=2)  # Sum over y dimension
    h_x = fit(Histogram, x, bins_x)
    px_individual = h_x.weights / (m * d_bins_x)
    
    # Marginal probabilities for Y
    py_marginal = sum(pxy, dims=1)  # Sum over x dimension
    h_y = fit(Histogram, y, bins_y)
    py_individual = h_y.weights / (m * d_bins_y)
    
    return pxy, px_individual, px_marginal, py_individual, py_marginal
end

"""
    _compute_mutual_information_integrand(
        pxy::Matrix,
        px_individual::Vector,
        py_individual::Vector,
        normalize_by_sample_size::Bool=false,
        sample_size::Int=0
    )::Matrix

Internal helper: compute the integrand for MI calculation.

# Arguments
- `pxy::Matrix`: Joint probability distribution
- `px_individual::Vector`: Marginal of X
- `py_individual::Vector`: Marginal of Y
- `normalize_by_sample_size::Bool`: Whether to normalize by sample size (for MIn)
- `sample_size::Int`: Sample size if normalizing

# Returns
- Matrix containing p(x,y) * log₂(p(x,y) / (p(x)p(y)))
"""
function _compute_mutual_information_integrand(
    pxy::Matrix,
    px_individual::Vector,
    py_individual::Vector,
    normalize_by_sample_size::Bool=false,
    sample_size::Int=0
)::Matrix
    
    # Initialize integrand array
    integrand = zeros(Float64, size(pxy))
    
    # Compute integrand value at each grid point
    for i in 1:length(px_individual)
        for j in 1:length(py_individual)
            # Skip if joint probability is zero (log(0) is undefined)
            if pxy[i, j] > 0
                # Compute log likelihood ratio
                log_ratio = log2(pxy[i, j] / (px_individual[i] * py_individual[j]))
                
                if normalize_by_sample_size
                    # MIn variant: include sample size normalization
                    integrand[i, j] = (pxy[i, j] / sample_size) * log_ratio
                else
                    # MI variant: standard integrand
                    integrand[i, j] = pxy[i, j] * log_ratio
                end
            end
        end
    end
    
    return integrand
end

end  # module MutualInformationAnalysis

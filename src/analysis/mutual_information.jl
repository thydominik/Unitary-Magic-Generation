"""
    MutualInformationAnalysis

Module for computing mutual information between two distributions using histogram-based estimation.

This module provides statistical analysis tools for quantifying dependence between variables
through information-theoretic measures. Two variants are provided:

1. **MI**: Fixed-range mutual information (assumes data in [0, 10])
2. **MIn**: Data-driven mutual information (normalizes to actual data range)

# Key Features

- Histogram-based probability density estimation from discrete samples
- 2D numerical integration using trapezoidal rule (from NumericalIntegration module)
- Proper normalization by bin widths for probability density
- Handles singular/zero probability cases to avoid log(0) errors
- Returns complete probability distributions for detailed analysis
- Type-stable implementation with proper type annotations

# Exported Functions

- `MI`: Compute mutual information with fixed range [0, 10]
- `MIn`: Compute mutual information with data-driven range

# Theory Background

Mutual Information (MI) quantifies the amount of information that one random variable 
carries about another, measured in bits.

## Mathematical Definition

```
MI(X;Y) = ∫∫ p(x,y) log₂(p(x,y) / (p(x)p(y))) dx dy
       = H(X) + H(Y) - H(X,Y)
```

where:
- p(x,y) is the joint probability density function
- p(x) and p(y) are the marginal probability densities
- H(X), H(Y), H(X,Y) are Shannon entropies (in bits)

## Interpretation

- MI = 0: X and Y are independent (no mutual information)
- MI > 0: X and Y are dependent (positive mutual information)
- MI increases with stronger statistical dependence
- Measured in bits when using log₂ (binary logarithm)
- MI is symmetric: MI(X;Y) = MI(Y;X)

# Implementation Details

## Histogram-Based Estimation

Since we only have discrete samples, we estimate the probability density using histograms:

1. Create 2D histogram of (x,y) samples partitioned into bins
2. Normalize histogram counts by (sample size × bin widths) to get probability densities
3. Compute marginals by summing over dimensions: p(x) = ∫ p(x,y) dy
4. Compute integrand = p(x,y) × log₂(p(x,y) / (p(x)p(y)))
5. Integrate using 2D trapezoidal rule

## Probability Density Normalization

The key to correct MI estimation is proper probability density normalization:

```
p_density(i,j) = histogram_count(i,j) / (N_samples × bin_width_x × bin_width_y)
```

This ensures the integral of p(x,y) over the entire domain equals 1:

```
∫∫ p(x,y) dx dy = 1
```

## Handling Zero Probabilities

The implementation carefully handles cases where p(x,y) = 0:
- Skips computation when p(x,y) ≤ 0 to avoid undefined log(0)
- Only computes log when probability is strictly positive
- This prevents NaN values in the integrand

# Example

```julia
using UnitaryMagicGeneration

# Create sample data with known dependence
x = randn(10000) .+ 5          # Normal distribution centered at 5
y = 0.7 .* x .+ 0.3 .* randn(10000)  # Y = 0.7*X + noise

# Compute mutual information with fixed range
mi, px_ind, px_marg, py_ind, py_marg, pxy, integrand = MI(x, y, 2^10)
println("Mutual Information: \$(mi) bits")
println("Positive MI indicates dependence between x and y")

# Compare with independent data
x_indep = randn(10000)
y_indep = randn(10000)
mi_indep, _ = MIn(x_indep, y_indep, 2^10)
println("MI for independent data: \$(mi_indep) (should be near 0)")
```
"""
module MutualInformationAnalysis

# ============================================================================
# MODULE IMPORTS
# ============================================================================
using StatsBase: fit, Histogram
using ..NumericalIntegration: double_integral_trapz

# ============================================================================
# EXPORTS
# ============================================================================
export MI, MIn

# ============================================================================
# PRIVATE HELPER FUNCTIONS
# ============================================================================

"""
    _compute_histogram_bins(
        data::AbstractVector{<:Real}, 
        n_bins::Int, 
        fixed_range::Union{Tuple{Real,Real}, Nothing}=nothing
    )::Tuple{AbstractRange, Float64}

Internal helper: compute histogram bins and bin width.

# Arguments
- `data::AbstractVector{<:Real}`: Data to bin
- `n_bins::Int`: Number of bins to create
- `fixed_range::Union{Tuple{Real,Real}, Nothing}`: Optional fixed range (min, max)

# Returns
- Tuple of (bins, bin_width) where bins is an AbstractRange and bin_width is Float64
"""
function _compute_histogram_bins(
    data::AbstractVector{<:Real}, 
    n_bins::Int, 
    fixed_range::Union{Tuple{Real,Real}, Nothing}=nothing
)::Tuple{AbstractRange, Float64}
    
    if fixed_range !== nothing
        # Use provided fixed range
        bins = range(fixed_range[1], fixed_range[2], n_bins)
    else
        # Use data-driven range from actual data
        bins = range(minimum(data), maximum(data), n_bins)
    end
    
    # Compute bin width (assumes uniform spacing)
    bin_width = Float64(bins[2] - bins[1])
    
    return bins, bin_width
end

"""
    _compute_probability_distributions(
        x::AbstractVector{<:Real},
        y::AbstractVector{<:Real},
        bins_x::AbstractRange,
        bins_y::AbstractRange
    )::Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}

Internal helper: compute joint and marginal probability distributions.

Computes normalized probability density functions from 2D histogram data.

# Returns
- Tuple of (pxy, px_individual, px_marginal, py_individual, py_marginal)
  - pxy: Joint probability distribution (normalized matrix)
  - px_individual: Individual x probabilities from histogram
  - px_marginal: Marginal probability p(x)
  - py_individual: Individual y probabilities from histogram
  - py_marginal: Marginal probability p(y)
"""
function _compute_probability_distributions(
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    bins_x::AbstractRange,
    bins_y::AbstractRange
)::Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}
    
    m = length(x)
    d_bins_x = Float64(bins_x[2] - bins_x[1])
    d_bins_y = Float64(bins_y[2] - bins_y[1])
    
    # Create 2D histogram and normalize to probability density
    joint_hist = fit(Histogram, (x, y), (bins_x, bins_y))
    pxy = joint_hist.weights ./ (m * d_bins_x * d_bins_y)
    
    # Marginal probabilities for X (sum over y dimension)
    px_marginal = vec(sum(pxy, dims=2))
    h_x = fit(Histogram, x, bins_x)
    px_individual = h_x.weights / (m * d_bins_x)
    
    # Marginal probabilities for Y (sum over x dimension)
    py_marginal = vec(sum(pxy, dims=1))
    h_y = fit(Histogram, y, bins_y)
    py_individual = h_y.weights / (m * d_bins_y)
    
    return pxy, px_individual, px_marginal, py_individual, py_marginal
end

"""
    _compute_mutual_information_integrand(
        pxy::Matrix{Float64},
        px_individual::Vector{Float64},
        py_individual::Vector{Float64},
        normalize_by_sample_size::Bool=false,
        sample_size::Int=0
    )::Matrix{Float64}

Internal helper: compute the integrand for MI calculation.

Computes p(x,y) * log₂(p(x,y) / (p(x)p(y))) at all grid points,
handling zero probabilities correctly.

# Arguments
- `pxy::Matrix{Float64}`: Joint probability distribution
- `px_individual::Vector{Float64}`: Marginal probability of X
- `py_individual::Vector{Float64}`: Marginal probability of Y
- `normalize_by_sample_size::Bool`: Whether to normalize by sample size (for MIn)
- `sample_size::Int`: Sample size if normalizing

# Returns
- Matrix{Float64}: The integrand p(x,y) * log₂(p(x,y) / (p(x)p(y)))
"""
function _compute_mutual_information_integrand(
    pxy::Matrix{Float64},
    px_individual::Vector{Float64},
    py_individual::Vector{Float64},
    normalize_by_sample_size::Bool=false,
    sample_size::Int=0
)::Matrix{Float64}
    
    # Initialize integrand array
    integrand = zeros(Float64, size(pxy))
    
    # Compute integrand value at each grid point
    for i in 1:length(px_individual)
        for j in 1:length(py_individual)
            # Skip if joint probability is zero (log(0) is undefined)
            if pxy[i, j] > 0.0
                # Compute log likelihood ratio
                px_py = px_individual[i] * py_individual[j]
                
                # Avoid division by zero
                if px_py > 0.0
                    log_ratio = log2(pxy[i, j] / px_py)
                    
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
    end
    
    return integrand
end

# ============================================================================
# PUBLIC FUNCTIONS
# ============================================================================

"""
    MI(
        x::AbstractVector{<:Real}, 
        y::AbstractVector{<:Real}, 
        N_bins::Int
    )::Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Matrix{Float64}, Matrix{Float64}}

Compute mutual information using fixed range [0, 10] with 2D trapezoidal integration.

This function assumes data falls within the range [0, 10]. For data outside this range,
use `MIn()` instead. This variant is useful when analyzing bounded quantities like 
probabilities or normalized statistics.

## Mathematical Background

Mutual Information is estimated as:

```
MI(X;Y) = ∫₀¹⁰ ∫₀¹⁰ p(x,y) log₂(p(x,y) / (p(x)p(y))) dx dy
```

where:
- p(x,y) is estimated from 2D histogram of samples
- Marginals p(x) and p(y) are obtained by summing over dimensions
- Integration performed using 2D trapezoidal rule

## Arguments

- `x::AbstractVector{<:Real}`: First data array (sample values)
  - Should contain values roughly in [0, 10]
  - Length can be arbitrary (typically 1000+ for stable estimation)
  - Type: any real numeric vector

- `y::AbstractVector{<:Real}`: Second data array (sample values)
  - Should contain values roughly in [0, 10]
  - Must have same length as x
  - Type: any real numeric vector

- `N_bins::Int`: Number of histogram bins (per dimension)
  - Controls resolution of probability density estimation
  - Typical range: 2¸ to 2¹⁴ (256 to 16384)
  - Larger = finer resolution but may overfit to noise
  - Rule of thumb: N_bins ≈ √(length(x))

## Returns

Tuple with 7 elements:

```julia
mi, px_individual, px_marginal, py_individual, py_marginal, pxy, integrand = MI(x, y, N_bins)
```

- `mi::Float64`: Mutual information value (in bits)
  - 0 for independent distributions
  - Positive for dependent distributions
  - Grows with stronger dependence

- `px_individual::Vector{Float64}`: Histogram counts for x (unnormalized)
  - Obtained by summing pxy over y dimension
  - Useful for analyzing marginal x distribution

- `px_marginal::Vector{Float64}`: Marginal probabilities for x
  - Normalized version of px_individual
  - Sums to 1 (integrated over domain)

- `py_individual::Vector{Float64}`: Histogram counts for y (unnormalized)
  - Obtained by summing pxy over x dimension

- `py_marginal::Vector{Float64}`: Marginal probabilities for y
  - Normalized version of py_individual
  - Sums to 1 (integrated over domain)

- `pxy::Matrix{Float64}`: Joint probability distribution p(x,y)
  - Size: (N_bins, N_bins)
  - Sums to 1 (integrated over entire domain)
  - Useful for analyzing dependence structure

- `integrand::Matrix{Float64}`: The integrand used in MI computation
  - Element-wise: p(x,y) * log₂(p(x,y) / (p(x)p(y)))
  - Useful for debugging and visualization of dependence regions

## Algorithm

1. Define bin edges: `range(0, 10, N_bins)`
2. Create 2D histogram of (x, y) samples
3. Normalize by sample size and bin widths to get probability densities
4. Compute marginals by dimension-wise summation
5. Compute integrand avoiding log(0) by checking `pxy[i,j] > 0`
6. Integrate using 2D trapezoidal rule (calls `double_integral_trapz`)
7. Return MI and all distributions for further analysis

## Example

```julia
using UnitaryMagicGeneration

# Create sample data with known dependence
x = randn(10000) .+ 5                      # Normal distribution centered at 5
y = 0.7 .* x .+ 0.3 .* randn(10000)       # Y correlated with X

# Compute mutual information with fixed range
mi, px_ind, px_marg, py_ind, py_marg, pxy, integrand = MI(x, y, 2^10)

println("Mutual Information: \$(mi) bits")
println("Higher MI indicates stronger dependence between x and y")

# Check with independent data
x_indep = randn(10000)
y_indep = randn(10000)
mi_indep, = MI(x_indep, y_indep, 2^10)
println("MI for independent data: \$(mi_indep) (should be near 0)")
```

## Notes

- **Fixed range assumption**: Assumes data in [0, 10]. Data outside will be truncated by histogram binning
- **Bin count selection**: Too few bins → information loss. Too many bins → overfitting to noise
- **Sample size**: More samples → more stable MI estimate
- **Smooth data**: Works best for smooth distributions; oscillatory data may need more bins
- **Zero handling**: Correctly handles p(x,y)=0 by skipping log(0) computation
"""
function MI(
    x::AbstractVector{<:Real}, 
    y::AbstractVector{<:Real}, 
    N_bins::Int
)::Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Matrix{Float64}, Matrix{Float64}}
    
    # ========================================================================
    # INPUT VALIDATION
    # ========================================================================
    if length(x) != length(y)
        error(
            "Length mismatch: x has length $(length(x)) but y has length $(length(y))"
        )
    end
    
    if N_bins < 2
        error("N_bins must be at least 2, got $N_bins")
    end
    
    # ========================================================================
    # STEP 1: CREATE HISTOGRAM BINS
    # ========================================================================
    # Define fixed bins in range [0, 10]
    bins_x = range(0, 10, N_bins)
    bins_y = range(0, 10, N_bins)
    
    # ========================================================================
    # STEP 2: COMPUTE PROBABILITY DISTRIBUTIONS
    # ========================================================================
    pxy, px_individual, px_marginal, py_individual, py_marginal = 
        _compute_probability_distributions(x, y, bins_x, bins_y)
    
    # ========================================================================
    # STEP 3: COMPUTE INTEGRAND
    # ========================================================================
    integrand = _compute_mutual_information_integrand(
        pxy, 
        px_individual, 
        py_individual,
        false,  # Don't normalize by sample size for MI
        0       # Sample size not needed for MI
    )
    
    # ========================================================================
    # STEP 4: COMPUTE 2D INTEGRAL USING TRAPEZOIDAL RULE
    # ========================================================================
    mi = double_integral_trapz(integrand, bins_x[1:end-1], bins_y[1:end-1])
    
    # ========================================================================
    # STEP 5: RETURN RESULTS
    # ========================================================================
    return mi, px_individual, px_marginal, py_individual, py_marginal, pxy, integrand
end

"""
    MIn(
        x::AbstractVector{<:Real}, 
        y::AbstractVector{<:Real}, 
        N_bins::Int
    )::Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Matrix{Float64}, Matrix{Float64}}

Compute mutual information using data-driven range (normalized MI).

Unlike `MI()`, this function determines bin ranges automatically from the data itself.
This is more robust for data with arbitrary ranges and non-uniform support.

## Mathematical Background

Normalized Mutual Information:

```
MI_n(X;Y) = ∫∫ p̃(x,y) log₂(p̃(x,y) / (p̃(x)p̃(y))) dx dy
```

where p̃ includes sample-size normalization factor.

## Differences from MI()

| Aspect | MI() | MIn() |
|--------|------|-------|
| **Range** | Fixed [0, 10] | Data-driven [min, max] |
| **Use case** | Bounded data | Arbitrary range |
| **Normalization** | Standard | Includes sample size factor |
| **Sample-size dependence** | Low | Moderate |

## Arguments

- `x::AbstractVector{<:Real}`: First data array (arbitrary range)
- `y::AbstractVector{<:Real}`: Second data array (arbitrary range)
- `N_bins::Int`: Number of histogram bins (see `MI()` for guidance)

## Returns

Same as `MI()`: tuple of (mi, px_individual, px_marginal, py_individual, py_marginal, pxy, integrand)

See `MI()` documentation for return value descriptions.

## Example

```julia
using UnitaryMagicGeneration

# Data with arbitrary range (not in [0, 10])
x = exp.(randn(5000))        # Exponential-like distribution
y = x .+ randn(5000)         # Dependent on x

# Use MIn for data-driven range
mi_norm, = MIn(x, y, 2^12)
println("Normalized MI: \$(mi_norm)")

# Compare with MI (if forcing [0,10] range)
mi_fixed, = MI(x, y, 2^12)  # Will truncate/distort data
println("Fixed-range MI: \$(mi_fixed)")
```

## Notes

- Recommended for real-world data with unknown or variable range
- More robust to outliers (only affects bin edges)
- May be slightly slower than `MI()` for very large datasets due to min/max computation
"""
function MIn(
    x::AbstractVector{<:Real}, 
    y::AbstractVector{<:Real}, 
    N_bins::Int
)::Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Matrix{Float64}, Matrix{Float64}}
    
    # ========================================================================
    # INPUT VALIDATION
    # ========================================================================
    if length(x) != length(y)
        error(
            "Length mismatch: x has length $(length(x)) but y has length $(length(y))"
        )
    end
    
    if N_bins < 2
        error("N_bins must be at least 2, got $N_bins")
    end
    
    # ========================================================================
    # STEP 1: CREATE DATA-DRIVEN HISTOGRAM BINS
    # ========================================================================
    # Determine bin ranges from actual data
    bins_x = range(minimum(x), maximum(x), N_bins)
    bins_y = range(minimum(y), maximum(y), N_bins)
    
    # ========================================================================
    # STEP 2: COMPUTE PROBABILITY DISTRIBUTIONS
    # ========================================================================
    pxy, px_individual, px_marginal, py_individual, py_marginal = 
        _compute_probability_distributions(x, y, bins_x, bins_y)
    
    # ========================================================================
    # STEP 3: COMPUTE INTEGRAND WITH SAMPLE-SIZE NORMALIZATION
    # ========================================================================
    m = length(x)
    integrand = _compute_mutual_information_integrand(
        pxy, 
        px_individual, 
        py_individual,
        true,   # Enable sample size normalization for MIn
        m       # Pass sample size
    )
    
    # ========================================================================
    # STEP 4: COMPUTE 2D INTEGRAL USING TRAPEZOIDAL RULE
    # ========================================================================
    mi = double_integral_trapz(integrand, bins_x[1:end-1], bins_y[1:end-1])
    
    # ========================================================================
    # STEP 5: RETURN RESULTS
    # ========================================================================
    return mi, px_individual, px_marginal, py_individual, py_marginal, pxy, integrand
end

end  # module MutualInformationAnalysis

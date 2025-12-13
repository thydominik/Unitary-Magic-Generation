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

# ============================================================================
# PUBLIC FUNCTIONS
# ============================================================================

"""
    MI(x::Vector, y::Vector, N_bins::Int)::Tuple

Compute mutual information using fixed range [0, 10] with trapezoidal integration.

This function assumes data falls within the range [0, 10]. For data outside this range,
use MIn() instead. This is useful when analyzing bounded quantities like probabilities
or normalized statistics.

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

- `x::Vector`: First data array (sample values)
  - Should contain values roughly in [0, 10]
  - Length can be arbitrary (typically 1000+)
  - Type: any numeric vector

- `y::Vector`: Second data array (sample values)
  - Should contain values roughly in [0, 10]
  - Must have same length as x
  - Type: any numeric vector

- `N_bins::Int`: Number of histogram bins
  - Controls resolution of probability estimation
  - Typical range: 2^8 to 2^14 (256 to 16384)
  - Larger = finer resolution but may suffer overfitting
  - Rule of thumb: N_bins ≈ sqrt(length(x))

## Returns

Tuple with 7 elements:

```julia
mi, px_individual, px_marginal, py_individual, py_marginal, pxy, integrand = MI(x, y, N_bins)
```

- `mi::Float64`: Mutual information value (in bits)
  - 0 for independent distributions
  - Positive for dependent distributions
  - Grows with stronger dependence

- `px_individual::Vector`: Marginal histogram counts for x
  - Obtained by summing pxy over y
  - Useful for analyzing x distribution

- `px_marginal::Vector`: Marginal probabilities for x
  - Normalized version of px_individual
  - Sum equals 1

- `py_individual::Vector`: Marginal histogram counts for y
  - Obtained by summing pxy over x

- `py_marginal::Vector`: Marginal probabilities for y
  - Normalized version of py_individual
  - Sum equals 1

- `pxy::Matrix`: Joint probability distribution p(x,y)
  - Size: (N_bins, N_bins)
  - Sum equals 1 (integrated over domain)
  - Useful for analyzing dependence structure

- `integrand::Matrix`: The integrand used in MI computation
  - p(x,y) * log₂(p(x,y) / (p(x)p(y)))
  - Useful for debugging and visualization

## Algorithm

1. Define bin edges: range(0, 10, N_bins)
2. Create 2D histogram of (x, y) samples
3. Normalize by sample size and bin widths to get probability densities
4. Compute marginals by dimension-wise summation
5. Compute integrand avoiding log(0) by checking pxy[i,j] > 0
6. Integrate using 2D trapezoidal rule (calls double_integral_trapz)
7. Return MI and all distributions for further analysis

## Example

```julia
using UnitaryMagic

# Create sample data
x = randn(10000) .+ 5  # Normal distribution centered at 5
y = 0.7 .* x .+ 0.3 .* randn(10000)  # Y correlated with X

# Compute mutual information
mi, px_ind, px_marg, py_ind, py_marg, pxy, integrand = MI(x, y, 2^10)

println("Mutual Information: $(mi) bits")
println("Higher MI indicates stronger dependence between x and y")

# Check if data is really independent
x_indep = randn(10000)
y_indep = randn(10000)
mi_indep, = MI(x_indep, y_indep, 2^10)
println("MI for independent data: $(mi_indep) (should be near 0)")
```

## Notes

- **Fixed range assumption**: Assumes data in [0, 10]. Data outside this range will be truncated by histogram binning
- **Bin count selection**: Too few bins → information loss. Too many bins → overfitting to noise
- **Sample size**: More samples → more stable MI estimate
- **Smooth data**: Works best for smooth distributions; oscillatory data may need more bins
- **Zero handling**: Correctly handles p(x,y)=0 by skipping log(0)

"""
function MI(
    x::AbstractVector, 
    y::AbstractVector, 
    N_bins::Int
)::Tuple{Float64, Vector, Vector, Vector, Vector, Matrix, Matrix}
    
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
    # Calls internal helper function
    pxy, px_individual, px_marginal, py_individual, py_marginal = 
        _compute_probability_distributions(x, y, bins_x, bins_y)
    
    # ========================================================================
    # STEP 3: COMPUTE INTEGRAND
    # ========================================================================
    # Calls internal helper function
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
    # IMPORTANT: This calls the function FROM NumericalIntegration module
    # NOT a duplicate implementation
    mi = double_integral_trapz(integrand, bins_x[1:end-1], bins_y[1:end-1])
    
    # ========================================================================
    # STEP 5: RETURN RESULTS
    # ========================================================================
    return mi, px_individual, px_marginal, py_individual, py_marginal, pxy, integrand
end


"""
    MIn(x::Vector, y::Vector, N_bins::Int)::Tuple

Compute mutual information using data-driven range (normalized MI).

Unlike MI(), this function determines bin ranges automatically from the data itself.
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
| **Sample-size dep.** | Low | Moderate |

## Arguments

- `x::Vector`: First data array (arbitrary range)
- `y::Vector`: Second data array (arbitrary range)
- `N_bins::Int`: Number of histogram bins (see MI() for guidance)

## Returns

Same as MI(): tuple of (mi, px_individual, px_marginal, py_individual, py_marginal, pxy, integrand)

## Example

```julia
using UnitaryMagic

# Data with arbitrary range (not in [0, 10])
x = exp.(randn(5000))  # Exponential-like distribution
y = x .+ randn(5000)   # Dependent on x

# Use MIn for data-driven range
mi_norm, = MIn(x, y, 2^12)
println("Normalized MI: $(mi_norm)")

# Compare with MI (if forcing [0,10] range)
mi_fixed, = MI(x, y, 2^12)  # Will truncate/distort data
println("Fixed-range MI: $(mi_fixed)")
```

## Notes

- Recommended for real-world data with unknown or variable range
- More robust to outliers (only affects bin edges)
- May be slower than MI() for very large datasets due to min/max computation
"""
function MIn(
    x::AbstractVector, 
    y::AbstractVector, 
    N_bins::Int
)::Tuple{Float64, Vector, Vector, Vector, Vector, Matrix, Matrix}
    
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
    # Calls internal helper function (same as MI)
    pxy, px_individual, px_marginal, py_individual, py_marginal = 
        _compute_probability_distributions(x, y, bins_x, bins_y)
    
    # ========================================================================
    # STEP 3: COMPUTE INTEGRAND WITH SAMPLE-SIZE NORMALIZATION
    # ========================================================================
    # Calls internal helper function with normalization enabled
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
    # IMPORTANT: This calls the function FROM NumericalIntegration module
    # NOT a duplicate implementation
    mi = double_integral_trapz(integrand, bins_x[1:end-1], bins_y[1:end-1])
    
    # ========================================================================
    # STEP 5: RETURN RESULTS
    # ========================================================================
    return mi, px_individual, px_marginal, py_individual, py_marginal, pxy, integrand
end

end  # module MutualInformationAnalysis

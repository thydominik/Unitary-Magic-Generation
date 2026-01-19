"""
    NumericalIntegration

Module for numerical integration utilities used throughout the package.

This module provides core numerical methods for computing integrals of discrete functions.
Primary use case: computing 2D integrals of probability distributions and information measures.

# Features
- Type-safe implementation with full error checking
- Clear, documented algorithms
- No external dependencies beyond Base Julia
- Suitable for smooth functions over regular grids
- Supports non-uniform grid spacing

# Exported Functions
- `double_integral_trapz`: Compute 2D definite integrals using trapezoidal rule

# Example

```julia
using UnitaryMagicGeneration.Analysis

# Define a test function: f(x,y) = exp(-(x² + y²))
x = range(0, 2, 100)
y = range(0, 2, 100)
f = [exp(-(x[i]^2 + y[j]^2)) for i in 1:100, j in 1:100]

# Compute the 2D integral
result = double_integral_trapz(f, x, y)
println("Integral ≈ \$result")
```
"""
module NumericalIntegration

# ============================================================================
# EXPORTS
# ============================================================================
# All functions that users should call are explicitly exported
export double_integral_trapz

# ============================================================================
# FUNCTION DEFINITIONS
# ============================================================================

"""
    double_integral_trapz(
        f::Array{<:Real,2}, 
        x::AbstractVector, 
        y::AbstractVector
    )::Float64

Computes a 2D double integral using the composite trapezoidal rule.

This function integrates a function f(x,y) over a 2D rectangular domain using
the trapezoidal rule for numerical integration. The function is evaluated on a
regular grid defined by x and y coordinate vectors.

## Algorithm

The 2D trapezoidal rule works as follows:
1. First, compute 1D trapezoidal integrals over x at each fixed y value
2. Then, compute a 1D trapezoidal integral of the results over y

This is equivalent to:
```
∫∫ f(x,y) dx dy ≈ Σᵢ Σⱼ (1/4) * dx * dy * [f(xᵢ,yⱼ) + f(xᵢ₊₁,yⱼ) + f(xᵢ,yⱼ₊₁) + f(xᵢ₊₁,yⱼ₊₁)]
```

## Arguments

- `f::Array{<:Real,2}`: 2D array containing function values at grid points
  - Dimensions: (nx, ny) where nx = length(x) and ny = length(y)
  - Typically computed as f[i,j] = func(x[i], y[j])
  - Must be real-valued (supports Int64, Float32, Float64, etc.)
  - Type stable: converts to Float64 internally

- `x::AbstractVector`: 1D vector of x-axis grid points
  - Length must equal first dimension of f
  - Should be sorted (ascending or descending, doesn't matter for trapz)
  - Can be Vector, UnitRange, LinRange, etc.
  - Grid spacing can be non-uniform

- `y::AbstractVector`: 1D vector of y-axis grid points
  - Length must equal second dimension of f
  - Should be sorted (ascending or descending)
  - Can be Vector, UnitRange, LinRange, etc.
  - Grid spacing can be non-uniform

## Returns

- `result::Float64`: The approximate value of the 2D integral
  - Always returns Float64 for type stability
  - Can be negative, zero, or positive
  - Accuracy depends on function smoothness and grid resolution

## Throws

- `DimensionMismatch`: If dimensions of f, x, and y don't match
  - Error message includes details about expected vs actual sizes
  - Helps identify whether problem is with f, x, or y

## Accuracy Considerations

- **Order of accuracy**: O(h²) where h is average grid spacing (second order)
- **Best for**: Smooth functions with regular or nearly-regular grid spacing
- **Performance vs accuracy**: Degrades for highly oscillatory functions
- **Grid density**: Generally need ~100+ points per axis for reasonable accuracy
- **Non-uniform grids**: Supported but may reduce accuracy slightly

## Implementation Details

The implementation uses the composite trapezoidal rule:

1. **X-Integration Phase**: For each y-value, compute ∫ f(x, y_j) dx
   - Uses formula: Σᵢ (1/2) * Δxᵢ * (f[i,j] + f[i+1,j])
   - Accumulates result in `integral_x[j]`

2. **Y-Integration Phase**: Integrate the x-integrated results over y
   - Uses formula: Σⱼ (1/2) * Δyⱼ * (integral_x[j] + integral_x[j+1])
   - Accumulates result in `integral_xy`

## Example

```julia
using UnitaryMagicGeneration

# Define a test function: f(x,y) = exp(-(x² + y²))
x = range(0, 2, 100)
y = range(0, 2, 100)
f = [exp(-(x[i]^2 + y[j]^2)) for i in 1:100, j in 1:100]

# Compute the 2D integral
result = double_integral_trapz(f, x, y)
println("Integral ≈ \$result")

# Compare with theoretical value for Gaussian
theory = π * exp(-2) * (2 * erf(2)/sqrt(π))^2
println("Theoretical ≈ \$theory")
println("Error: \$(abs(result - theory))")
```

## Performance Notes

- Memory efficient: Uses O(ny) temporary storage for x-integration results
- No allocation of large intermediate 2D arrays
- Type-stable: always returns Float64
- Supports AbstractVector types (works with views efficiently)
"""
function double_integral_trapz(
    f::Array{<:Real,2}, 
    x::AbstractVector, 
    y::AbstractVector
)::Float64
    
    # ========================================================================
    # INPUT VALIDATION
    # ========================================================================
    # Check that dimensions match before proceeding
    nx, ny = size(f)
    
    if length(x) != nx || length(y) != ny
        error(
            "Dimension mismatch in double_integral_trapz:\n" *
            "  Function array f has shape: ($nx, $ny)\n" *
            "  X-axis vector has length: $(length(x))\n" *
            "  Y-axis vector has length: $(length(y))\n" *
            "  Expected: f with shape ($(length(x)), $(length(y)))"
        )
    end

    # ========================================================================
    # STEP 1: INTEGRATE OVER X FOR EACH FIXED Y
    # ========================================================================
    # For each y value, compute ∫ f(x, y) dx using 1D trapezoidal rule
    # This produces a 1D array of results, one for each y value
    integral_x = zeros(Float64, ny)  # Results of x-integration at each y
    
    for j in 1:ny
        # Integrate along x at this fixed y value
        for i in 1:(nx-1)
            # Trapezoidal rule: (1/2) * Δx * (f₀ + f₁)
            dx = Float64(x[i+1] - x[i])
            f_i_j = Float64(f[i, j])
            f_ip1_j = Float64(f[i+1, j])
            integral_x[j] += 0.5 * dx * (f_i_j + f_ip1_j)
        end
    end

    # ========================================================================
    # STEP 2: INTEGRATE OVER Y OF THE PREVIOUS RESULTS
    # ========================================================================
    # Now integrate the x-integrated values over y
    integral_xy = 0.0
    
    for j in 1:(ny-1)
        # Trapezoidal rule: (1/2) * Δy * (I₀ + I₁)
        dy = Float64(y[j+1] - y[j])
        integral_xy += 0.5 * dy * (integral_x[j] + integral_x[j+1])
    end

    # ========================================================================
    # RETURN RESULT
    # ========================================================================
    return integral_xy
end

end  # module NumericalIntegration

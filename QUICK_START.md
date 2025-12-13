# Quick Start Guide: Refactored UnitaryMagic Package

## Overview

This guide will get you up and running with the refactored modular structure in under 5 minutes.

## Installation

### Option 1: Using the Modules Directly

```bash
cd Unitary-Magic-Generation
cd refactor/modular-structure  # or switch to branch
```

```julia
# Add to your Julia script
include("src/UnitaryMagic.jl")
using .UnitaryMagic
```

### Option 2: As a Julia Package (Recommended for Future)

```bash
julia> ]
pkg> develop /path/to/Unitary-Magic-Generation

# In any Julia session:
julia> using UnitaryMagic
```

---

## Basic Usage

### Computing Mutual Information

```julia
using UnitaryMagic
using Random

# Create sample data
Random.seed!(42)
x = randn(1000) .+ 5  # Magic distribution
y = randn(1000) .+ 5  # Entanglement distribution

# Compute mutual information (fixed range [0, 10])
mi, px_ind, px_marg, py_ind, py_marg, pxy, integrand = MI(x, y, 2^10)

println("Mutual Information: $mi")
println("P(x) shape: $(size(px_ind))")
println("P(x,y) shape: $(size(pxy))")
```

### Using Normalized MI

```julia
# Data-driven range (normalized)
mi_norm, = MIn(x, y, 2^10)

println("Normalized MI: $mi_norm")
```

### Using the 2D Integration Utility

```julia
using UnitaryMagic

# Create test function
f = exp.(-randn(100, 100).^2)
x = range(0, 1, 100)
y = range(0, 1, 100)

# Compute 2D integral
result = double_integral_trapz(f, x, y)
println("Integral result: $result")
```

---

## Common Workflows

### Workflow 1: Single MI Calculation

```julia
using UnitaryMagic, JLD2

# Load your data
data = load("my_data.jld2")
x = data["magic_values"]
y = data["entanglement_values"]

# Compute MI with 4096 bins
mi, = MI(x, y, 2^12)

println("MI = $mi bits")
```

### Workflow 2: Sample Size Dependence Study

```julia
using UnitaryMagic, Random, ProgressBars

x = randn(10000)
y = randn(10000)

MI_values = Float64[]
Sample_sizes = Int[]

for m in ProgressBar(0:14)
    sample_size = 2^m
    if sample_size ≤ length(x)
        indices = randperm(length(x))[1:sample_size]
        mi, = MI(x[indices], y[indices], 2^12)
        push!(MI_values, mi)
        push!(Sample_sizes, sample_size)
    end
end

println("Sample sizes: $Sample_sizes")
println("MI values: $MI_values")
```

### Workflow 3: Batch Analysis

```julia
include("examples/mutual_information_analysis.jl")

# Use provided functions
results = Dict()
for n in 2:10
    x = load("data_N_$n.jld2")["magic"]
    y = load("data_N_$n.jld2")["entanglement"]
    
    mi_vals, sample_sizes = analyze_sample_size_dependence(x, y)
    results[n] = (mi_vals, sample_sizes)
end

# Plot results
plot_sample_size_scaling(results, 2:10)
plot_mi_vs_qubit_count(results, 20)
```

---

## API Reference

### Core Functions

#### `MI(x, y, N_bins)` - Fixed-Range Mutual Information

```julia
mi, px_ind, px_marg, py_ind, py_marg, pxy, integrand = MI(x, y, N_bins)
```

**Arguments**:
- `x::Vector`: First data array
- `y::Vector`: Second data array
- `N_bins::Int`: Number of histogram bins

**Returns**:
- `mi::Float64`: Mutual information value (bits)
- `px_ind, py_ind`: Individual marginal distributions
- `px_marg, py_marg`: Marginal probabilities
- `pxy::Matrix`: Joint probability distribution
- `integrand::Matrix`: The integrand used in computation

**Notes**:
- Uses fixed range [0, 10] for both axes
- Returns 0 for independent distributions
- Typically use `2^10` to `2^14` bins

#### `MIn(x, y, N_bins)` - Data-Driven MI

```julia
mi_norm, px_ind, px_marg, py_ind, py_marg, pxy, integrand = MIn(x, y, N_bins)
```

**Arguments**: Same as `MI()`

**Returns**: Same as `MI()`

**Notes**:
- Uses data range [min(x), max(x)] and [min(y), max(y)]
- Better for data with non-uniform support
- More sample-size dependent

#### `double_integral_trapz(f, x, y)` - 2D Trapezoidal Integration

```julia
result = double_integral_trapz(f, x, y)
```

**Arguments**:
- `f::Matrix`: 2D array of function values
- `x::Vector`: x-axis grid points
- `y::Vector`: y-axis grid points

**Returns**:
- `result::Float64`: Computed integral value

**Notes**:
- `size(f)` must match `length(x)` × `length(y)`
- Throws error if dimensions don't match

---

## Module Structure

```
├── src/
│   ├── UnitaryMagic.jl                    # Main module (include this)
│   ├── utils/
│   │   └── numerical_integration.jl       # Integration utilities
│   └── analytics/
│       └── mutual_information.jl          # MI computation
├── examples/
│   └── mutual_information_analysis.jl     # Analysis pipeline
├── REFACTORING_GUIDE.md                   # Detailed architecture
├── CODE_REVIEW.md                         # Code review & comments
└── QUICK_START.md                         # This file
```

**To use everything**:
```julia
include("src/UnitaryMagic.jl")
using .UnitaryMagic
# All functions available: MI, MIn, double_integral_trapz
```

---

## Examples

### Example 1: Independent vs Dependent Variables

```julia
using UnitaryMagic
using Plots

# Create data
x = randn(5000)
y_indep = randn(5000)      # Independent
y_dep = x .+ 0.1 .* randn(5000)  # Dependent

mi_indep, = MI(x, y_indep, 1000)
mi_dep, = MI(x, y_dep, 1000)

println("MI(independent) = $mi_indep")
println("MI(dependent) = $mi_dep")

# Plot
scatter(x, y_indep, label="Independent")
scatter!(x, y_dep, label="Dependent")
```

### Example 2: Gaussian Distribution Validation

```julia
using UnitaryMagic, Distributions

# Create correlated Gaussians
μ = 5.0
σ = 0.5
n = 10000

ρ = 0.7  # Correlation coefficient
x = randn(n) .* σ .+ μ
y = ρ .* x .+ sqrt(1 - ρ^2) .* randn(n) .* σ .+ μ

mi_computed, = MI(x, y, 2^12)

# Theoretical: MI = -0.5 * log(1 - ρ²)
mi_theory = -0.5 * log2(1 - ρ^2)

println("Computed MI: $mi_computed")
println("Theoretical MI: $mi_theory")
println("Error: $(abs(mi_computed - mi_theory))")
```

### Example 3: Convergence with Sample Size

```julia
using UnitaryMagic, Plots, Random

Random.seed!(42)
n_total = 100000
x = randn(n_total)
y = x .+ 0.1 .* randn(n_total)

mi_vals = Float64[]
sizes = Int[]

for log_m in 0:16
    m = min(2^log_m, n_total)
    indices = randperm(n_total)[1:m]
    mi, = MI(x[indices], y[indices], 2^10)
    push!(mi_vals, mi)
    push!(sizes, m)
end

plot(log2.(sizes), mi_vals,
    xlabel="Sample Size - log₂(N)",
    ylabel="Mutual Information (bits)",
    legend=false,
    title="MI Convergence",
    markerstrokewidth=0,
    markersize=6)
```

---

## Performance Tips

### 1. Choose Bin Count Wisely

```julia
# Too few bins: information loss
mi_coarse, = MI(x, y, 100)

# Good balance: computational cost vs accuracy
mi_good, = MI(x, y, 2^12)  # 4096 bins

# Too many bins: overfitting to sample noise
mi_fine, = MI(x, y, 2^16)  # 65536 bins (slow)

# Rule of thumb: bins ~ O(sqrt(N_samples))
N_bins_suggested = ceil(Int, sqrt(length(x)))
```

### 2. Cache Results for Repeated Use

```julia
using JLD2

# Compute once
mi, px_ind, px_marg, py_ind, py_marg, pxy, integrand = MI(x, y, 2^12)

# Save
@save "mi_cache.jld2" mi px_ind px_marg py_ind py_marg pxy integrand

# Later: load
@load "mi_cache.jld2" mi px_ind px_marg py_ind py_marg pxy integrand
```

### 3. Use Views for Large Arrays

```julia
using UnitaryMagic

large_x = randn(1_000_000)
large_y = randn(1_000_000)

# Create view (no copy)
slice_x = @view large_x[1:10000]
slice_y = @view large_y[1:10000]

mi, = MI(slice_x, slice_y, 2^10)  # Much faster
```

---

## Troubleshooting

### Issue: "Type error in function application"

**Cause**: Wrong argument types
```julia
# WRONG
mi = MI([1, 2, 3], [4, 5, 6], "100")  # N_bins should be Int

# CORRECT
mi, = MI([1.0, 2.0, 3.0], [4.0, 5.0, 6.0], 100)
```

### Issue: "Size of f must match lengths..."

**Cause**: Dimension mismatch in integration
```julia
# WRONG
f = randn(10, 20)
x = range(0, 1, 10)
y = range(0, 1, 10)  # Should be 20
double_integral_trapz(f, x, y)

# CORRECT
f = randn(10, 20)
x = range(0, 1, 10)
y = range(0, 1, 20)
double_integral_trapz(f, x, y)
```

### Issue: Very small or zero MI values

**Possible causes**:
1. Data is actually independent
2. Bin count too high (overfitting)
3. Data ranges don't fit [0, 10] for fixed MI

**Solutions**:
```julia
# Check data ranges
print("x range: [$(minimum(x)), $(maximum(x))]")
print("y range: [$(minimum(y)), $(maximum(y))]")

# Use normalized MI if data outside [0, 10]
mi_norm, = MIn(x, y, 2^12)

# Try fewer bins
mi_coarse, = MI(x, y, 2^8)
```

---

## Next Steps

1. ✅ Read [REFACTORING_GUIDE.md](REFACTORING_GUIDE.md) for architecture details
2. 🔲 See [CODE_REVIEW.md](CODE_REVIEW.md) for detailed comments
3. 🔲 Run `examples/mutual_information_analysis.jl` for full pipeline
4. 🔲 Check out existing examples in `Examples/` directory

---

## Getting Help

- **Module Documentation**: Run `? MI` in Julia REPL
- **Examples**: See `examples/` directory
- **Architecture**: Read REFACTORING_GUIDE.md
- **Code Comments**: Check CODE_REVIEW.md

---

**Version**: Refactor Branch  
**Last Updated**: 2025-12-13

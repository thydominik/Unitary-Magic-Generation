"""
    UnitaryMagic

A comprehensive Julia package for analyzing magic content and entanglement in unitary circuits.

This is the main entry point for the UnitaryMagic package. All public API functions are
accessed through this module.

# Core Modules

## NumericalIntegration
Provides numerical integration utilities:
- `double_integral_trapz`: 2D integration using trapezoidal rule

## MutualInformationAnalysis  
Provides information-theoretic analysis:
- `MI`: Compute mutual information with fixed range [0, 10]
- `MIn`: Compute mutual information with data-driven range

# Features

- Type-safe function signatures with full documentation
- Comprehensive error handling with informative messages
- No global state - pure functions only
- Fully modular design - each component independent

# Usage Examples

## Basic MI Computation
```julia
using UnitaryMagic

x = randn(1000)
y = randn(1000)
mi, px_ind, px_marg, py_ind, py_marg, pxy, integrand = MI(x, y, 2^10)
println("Mutual Information: \$mi bits")
```

## Using Data-Driven MI
```julia
mi_norm, = MIn(x, y, 2^12)
```

## 2D Integration
```julia
f = randn(100, 100)
x_grid = range(0, 1, 100)
y_grid = range(0, 1, 100)
result = double_integral_trapz(f, x_grid, y_grid)
```

# Package Structure

```
UnitaryMagic/
├─ UnitaryMagic.jl (this file - main entry point)
├─ utils/
│  └─ numerical_integration.jl (integration utilities)
└─ analytics/
   └─ mutual_information.jl (MI computation)
```
"""
module UnitaryMagic

# ============================================================================
# MODULE STRUCTURE AND INCLUDES
# ============================================================================
# This module follows a hierarchical structure where:
# 1. Utility modules are included first (they have no dependencies)
# 2. Analytics modules are included second (they depend on utilities)
# 3. All public functions are re-exported for user convenience

# Include utility modules (base functionality)
include("utils/numerical_integration.jl")
using .NumericalIntegration

# Include analytics modules (depends on utilities)
include("analytics/mutual_information.jl")
using .MutualInformationAnalysis

# ============================================================================
# PUBLIC API EXPORTS
# ============================================================================
# All exported functions are called FROM their respective modules
# Users should access all functions through the UnitaryMagic namespace

export 
    # From NumericalIntegration module
    double_integral_trapz,
    # From MutualInformationAnalysis module
    MI, 
    MIn

end  # module UnitaryMagic

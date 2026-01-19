"""
    Utils

Utility functions for UnitaryMagicGeneration package.

This module provides general-purpose helper functions and utilities used throughout
the package, including numerical integration tools.

# Submodules

- `NumericalIntegration`: Numerical integration utilities

# Exported Functions

- `double_integral_trapz`: Compute 2D definite integrals using trapezoidal rule

# Example

```julia
using UnitaryMagicGeneration.Utils

# Define a test function: f(x,y) = exp(-(x² + y²))
x = range(0, 2, 100)
y = range(0, 2, 100)
f = [exp(-(x[i]^2 + y[j]^2)) for i in 1:100, j in 1:100]

# Compute the 2D integral
result = double_integral_trapz(f, x, y)
println("Integral ≈ \$result")
```
"""
module Utils

# ============================================================================
# SUBMODULE INCLUDES
# ============================================================================
include("numerical_integration.jl")

# ============================================================================
# EXPORTS
# ============================================================================
# Re-export all public functions from submodules
export double_integral_trapz

# Also export the submodule itself for direct access
export NumericalIntegration

end  # module Utils

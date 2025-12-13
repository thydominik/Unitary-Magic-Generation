"""
    UnitaryMagic

A comprehensive Julia package for analyzing magic content and entanglement in unitary circuits.

This module provides tools for:
- Computing and analyzing magic states in quantum systems
- Measuring entanglement properties
- Statistical analysis via mutual information
- Random unitary generation and sampling

# Modules
- `NumericalIntegration`: Numerical integration utilities
- `MutualInformationAnalysis`: Mutual information computation and analysis

# Usage
```julia
using UnitaryMagic

# Example: Computing mutual information
x = randn(1000)
y = randn(1000)
mi, px_ind, px_marg, py_ind, py_marg, pxy, integrand = MI(x, y, 100)
```
"""
module UnitaryMagic

include("utils/numerical_integration.jl")
using .NumericalIntegration

include("analytics/mutual_information.jl")
using .MutualInformationAnalysis

# Re-export public functions
export double_integral_trapz, MI, MIn

end  # module UnitaryMagic

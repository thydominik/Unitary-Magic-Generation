"""
    Analysis

Analysis and post-processing functionality for UnitaryMagicGeneration package.

This module provides tools for statistical analysis of quantum circuits and states,
including mutual information calculations and other information-theoretic measures.

# Submodules

- `MutualInformationAnalysis`: Compute mutual information between distributions

# Exported Functions

All functions from submodules are automatically available through Analysis module.
See individual submodule documentation for details.

# Example

```julia
using UnitaryMagicGeneration.Analysis

# Create correlated data
x = randn(10000) .+ 5
y = 0.7 .* x .+ 0.3 .* randn(10000)

# Compute mutual information using the MutualInformationAnalysis submodule
mi, = MutualInformationAnalysis.MI(x, y, 2^10)
println("Mutual Information: \$mi bits")
```
"""
module Analysis

# ============================================================================
# SUBMODULE INCLUDES
# ============================================================================
# Import NumericalIntegration for use by MutualInformationAnalysis
using ..NumericalIntegration

# Include analysis submodules
include("mutual_information.jl")

# ============================================================================
# EXPORTS
# ============================================================================
# Re-export all public functions from submodules
export MI, MIn

# Also export the submodule itself for direct access
export MutualInformationAnalysis

end  # module Analysis

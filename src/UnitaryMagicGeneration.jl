"""
    UnitaryMagicGeneration

A Julia package for studying magic state generation and unitary circuits in quantum computing.

This package provides tools for analyzing and generating unitary circuits with specific magic 
state properties, including dissipative dynamics, stabilizer code operations, and resource
theory applications.

# Main Modules

- `Analysis`: Statistical analysis and information-theoretic measures
  - Mutual information computation
  - Probability distribution analysis
- `Utils`: Utility functions for numerical computation
  - Numerical integration tools
- `Core`: Core quantum computing algorithms (in development)
- `Circuits`: Circuit generation and analysis (in development)

# Quick Start

```julia
using UnitaryMagicGeneration

# Compute mutual information between two data distributions
x = randn(5000) .+ 5
y = 0.7 .* x .+ 0.3 .* randn(5000)

mi, = MI(x, y, 2^10)  # Fixed-range [0, 10] MI
println("Mutual Information: \$mi bits")

mi_norm, = MIn(x, y, 2^10)  # Data-driven MI
println("Normalized MI: \$mi_norm bits")
```

# Module Organization

The package is organized into functional submodules:

- **Analysis**: Information-theoretic measures and statistical analysis
  - MI(): Mutual information with fixed range
  - MIn(): Mutual information with data-driven range
  - Supporting probability distribution analysis

- **Utils**: General-purpose utilities
  - double_integral_trapz(): 2D numerical integration
  - Supporting numerical functions

- **Core**: Core quantum computing functionality (scaffold)
  - To be populated with core algorithms

- **Circuits**: Circuit-related functionality (scaffold)
  - To be populated with circuit generation and analysis

# Features

- Unitary circuit generation and analysis
- Magic state quantification and analysis
- Stabilizer circuit operations
- Dissipative system dynamics
- Clifford circuit analysis
- Information-theoretic measures

# Documentation

For detailed documentation, see the module-specific help:

```julia
?Analysis
?Analysis.MI
?Analysis.MIn
?Utils
?Utils.double_integral_trapz
```

# Examples

See the `examples/` directory for example scripts and notebooks.

"""
module UnitaryMagicGeneration

# ============================================================================
# MODULE INCLUDES
# ============================================================================
# Include submodules in dependency order
include("utils/Utils.jl")
include("analysis/Analysis.jl")
include("core/Core.jl")
include("circuits/Circuits.jl")

# ============================================================================
# EXPORTS
# ============================================================================
# Export all public API functions and submodules

# Analysis submodule and its functions
export Analysis, MI, MIn

# Utils submodule and its functions
export Utils, double_integral_trapz

# Core submodule (scaffold)
export Core

# Circuits submodule (scaffold)
export Circuits

end # module UnitaryMagicGeneration

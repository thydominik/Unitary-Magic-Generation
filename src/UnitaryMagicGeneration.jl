"""    UnitaryMagicGeneration

A Julia package for studying magic state generation and unitary circuits in quantum computing.

This package provides tools for analyzing and generating unitary circuits with specific magic 
state properties, including dissipative dynamics, stabilizer code operations, and resource
theory applications.

# Features

- Unitary circuit generation and analysis
- Magic state quantification
- Stabilizer circuit operations
- Dissipative system dynamics
- Clifford circuit analysis

# Quick Start

```julia
using UnitaryMagicGeneration

# Generate a random unitary circuit
# (specific functions depend on your core API)
```

See the documentation for detailed usage examples.
"""
module UnitaryMagicGeneration

export 
    # Core functions (to be populated after code refactoring)
    # Example placeholders:
    # generate_unitary,
    # calculate_magic,
    # analyze_circuit

include("core/Core.jl")
include("utils/Utils.jl")
include("analysis/Analysis.jl")
include("circuits/Circuits.jl")

end # module

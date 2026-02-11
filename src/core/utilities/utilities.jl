"""
    utilities

General helper utilities used across the repository.

This module is intentionally small and contains only general-purpose helper functions.

All identifiers and docstrings are ASCII-only.
"""
module utilities

include("numerical_integration.jl")

using .numerical_integration

export double_integral_trapz
export numerical_integration

end # module utilities

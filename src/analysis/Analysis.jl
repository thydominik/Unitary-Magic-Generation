"""
    analysis

Convenience module for analysis helpers.

This repository is organized as a project-style codebase, not a registered Julia package.
The analysis code in `src/analysis` is intended for scripts, notebooks, and one-off studies.

This module exists only as a thin convenience layer that makes it easy to access selected
core utilities from analysis code without duplicating implementations.

Notes
- All identifiers and docstrings are ASCII-only.
- This module intentionally reuses implementations from `src/core`.
"""
module analysis

# Core utilities commonly used from analysis.
include(joinpath(@__DIR__, "..", "core", "utilities", "utilities.jl"))
include(joinpath(@__DIR__, "..", "core", "mutual_information", "mutual_information.jl"))

using .utilities
using .mutual_information

export utilities,
    mutual_information

export double_integral_trapz
export mutual_information_discrete,
    mutual_information_histogram

end # module analysis

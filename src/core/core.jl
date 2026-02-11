"""
    core

Entry-point module that aggregates the repository's core building blocks.

This repository is not a registered Julia package, so the code is organized as a project-style
codebase under `src/core`. This file exists purely as a convenience module that:
- Includes the core submodules in a single place.
- Re-exports the most commonly used functions.

Notes
- All identifiers and docstrings are ASCII-only.
- The canonical implementations live in the snake_case files within each subfolder.
- Deprecated wrapper files may exist to preserve old include paths while refactoring.
"""
module core

# Canonical submodules (snake_case)
include(joinpath(@__DIR__, "magic", "magic.jl"))
include(joinpath(@__DIR__, "random_unitaries", "random_unitaries.jl"))
include(joinpath(@__DIR__, "circuits", "circuits.jl"))
include(joinpath(@__DIR__, "entanglement", "entanglement.jl"))
include(joinpath(@__DIR__, "utilities", "utilities.jl"))
include(joinpath(@__DIR__, "mutual_information", "mutual_information.jl"))

using .magic
using .random_unitaries
using .circuits
using .entanglement
using .utilities
using .mutual_information

export magic,
    random_unitaries,
    circuits,
    entanglement,
    utilities,
    mutual_information

# Re-export selected functions (keep this list short and stable)
export generate_all_pauli_strings,
    pauli_matrix,
    pauli_operator_list,
    measure_magic_pure,
    measure_magic_mixed

export cue_matrix,
    random_unitary_matrix

export brickwall_layer_matrix,
    brickwall_unitary_matrix

export reduced_density_matrix,
    von_neumann_entropy,
    renyi_entropy,
    calculate_entanglement,
    calculate_renyi_entropy,
    entanglement_negativity,
    logarithmic_negativity

export double_integral_trapz

export mutual_information_discrete,
    mutual_information_histogram

end # module core

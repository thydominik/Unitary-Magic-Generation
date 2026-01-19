"""
    Core

Core algorithms and quantum primitives for UnitaryMagicGeneration.

This module provides the main algorithmic components of the package, including:

- Magic quantification routines (`Magic` submodule)
- Random unitary and brick-wall circuit generation (`RandomUnitaries` submodule)

Additional core functionality will be added here as the refactoring progresses.

# Submodules

- `Magic`: Pauli operators and magic quantification
- `RandomUnitaries`: Haar-random unitaries and brick-wall circuits

# Examples

```julia
using UnitaryMagicGeneration

# Generate a random brick-wall circuit on 4 qubits
U = generate_brickwall_unitary_circuit(4)

# Construct Pauli operators and measure magic of a random state
using LinearAlgebra

n = 2
ψ = randn(ComplexF64, 2^n)
ψ /= norm(ψ)

pauli_strings = generate_all_pauli_strings(n)
pauli_ops = pauli_operator_list(pauli_strings, n)
magic, Ξ = measure_magic_pure(ψ, pauli_ops)
```
"""
module Core

# Submodules
include("Magic.jl")
include("RandomUnitaries.jl")

# Re-export key functions from submodules
export 
    Magic,
    RandomUnitaries,
    pauli_operator_list,
    pauli_matrix,
    measure_magic_pure,
    measure_magic_mixed,
    generate_all_pauli_strings,
    cue_matrix,
    generate_regular_unitary_circuit,
    generate_brickwall_unitary_circuit

end # module Core

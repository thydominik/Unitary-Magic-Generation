"""
    Core

Core algorithms and quantum primitives for UnitaryMagicGeneration.

This module provides the main algorithmic components of the package, including:

- Magic quantification routines (`Magic` submodule)
- Random unitary and brick-wall circuit generation (`RandomUnitaries` submodule)
- Entanglement quantification (`Entanglement` submodule)
- Rényi entropy computation (`Entropy` submodule)

# Submodules

- `Magic`: Pauli operators and magic quantification
- `RandomUnitaries`: Haar-random unitaries and brick-wall circuits
- `Entanglement`: Reduced density matrices and entanglement entropy
- `Entropy`: Rényi entropy across bipartitions

# Examples

```julia
using UnitaryMagicGeneration
using LinearAlgebra

# Generate a random brick-wall circuit on 4 qubits
U = generate_brickwall_unitary_circuit(4)

# Construct Pauli operators and measure magic of a random state
n = 2
ψ = randn(ComplexF64, 2^n)
ψ /= norm(ψ)

pauli_strings = generate_all_pauli_strings(n)
pauli_ops = pauli_operator_list(pauli_strings, n)
magic, Ξ = measure_magic_pure(ψ, pauli_ops)

# Compute entanglement entropy across bipartitions
entropies, subsystems = calculate_entanglement(ψ; subsystems="All")

# Compute Rényi-2 entropy (collision entropy)
S2 = calculate_renyi_entropy(ψ, 2.0; subsystems="Middle")
```
"""
module Core

# Submodules
include("Magic.jl")
include("RandomUnitaries.jl")
include("Entanglement.jl")
include("Entropy.jl")

# Re-export key functions from submodules
export 
    # Submodule names
    Magic,
    RandomUnitaries,
    Entanglement,
    Entropy,
    
    # Magic quantification
    pauli_operator_list,
    pauli_matrix,
    measure_magic_pure,
    measure_magic_mixed,
    generate_all_pauli_strings,
    
    # Random unitaries
    cue_matrix,
    generate_regular_unitary_circuit,
    generate_brickwall_unitary_circuit,
    
    # Entanglement
    reduced_density_matrix,
    von_neumann_entropy,
    calculate_entanglement,
    
    # Entropy
    renyi_entropy,
    calculate_renyi_entropy

end # module Core

"""
    RandomUnitaries

Core routines for generating random unitary matrices and brick-wall unitary circuits.

This module provides:
- Haar-random unitaries from the Circular Unitary Ensemble (CUE)
- Random brick-wall unitary circuits with alternating 2-qubit layers

# Exported API

- `cue_matrix` — Haar-random unitary from CUE
- `generate_regular_unitary_circuit` — single-layer CUE unitary on `n_qubits`
- `generate_brickwall_unitary_circuit` — depth-`D` brick-wall random unitary circuit
"""
module RandomUnitaries

using LinearAlgebra
using RandomMatrices

export cue_matrix, generate_regular_unitary_circuit, generate_brickwall_unitary_circuit

"""
    cue_matrix(dim::Integer; beta::Integer=2) -> Matrix{ComplexF64}

Draw a Haar-random unitary matrix from the Circular Unitary Ensemble (CUE).

# Arguments
- `dim`: Dimension of the Hilbert space (e.g. `2^n_qubits`)
- `beta`: Dyson index (default 2 for complex unitary ensemble)

# Returns
- `U::Matrix{ComplexF64}`: Haar-random unitary of size `dim × dim`
"""
function cue_matrix(dim::Integer; beta::Integer=2)::Matrix{ComplexF64}
    dim ≤ 0 && throw(ArgumentError("dim must be positive, got $dim"))
    U = rand(Haar(beta), dim)
    return Matrix{ComplexF64}(U)
end

"""
    generate_regular_unitary_circuit(n_qubits::Integer) -> Matrix{ComplexF64}

Generate a single-layer Haar-random unitary acting on `n_qubits`.

This is simply a CUE matrix on the full `2^n_qubits`-dimensional Hilbert space.
"""
function generate_regular_unitary_circuit(n_qubits::Integer)::Matrix{ComplexF64}
    n_qubits ≤ 0 && throw(ArgumentError("n_qubits must be positive, got $n_qubits"))
    dim = 2^n_qubits
    return cue_matrix(dim)
end

"""
    generate_brickwall_unitary_circuit(n_qubits::Integer; depth::Integer=5*n_qubits) -> Matrix{ComplexF64}

Generate a brick-wall random unitary circuit of given depth.

The circuit consists of alternating layers of 2-qubit unitaries sampled from CUE(4):
- Odd layers act on qubits (1,2), (3,4), ...
- Even layers act on qubits (2,3), (4,5), ...

Open boundary conditions are used; if a pair does not fit, an identity is applied.

# Arguments
- `n_qubits`: Number of qubits in the circuit
- `depth`: Number of layers (default `5*n_qubits`)

# Returns
- `U::Matrix{ComplexF64}`: Unitary matrix representing the full circuit
"""
function generate_brickwall_unitary_circuit(
    n_qubits::Integer;
    depth::Integer=5*n_qubits,
)::Matrix{ComplexF64}
    n_qubits ≤ 0 && throw(ArgumentError("n_qubits must be positive, got $n_qubits"))
    depth ≤ 0 && throw(ArgumentError("depth must be positive, got $depth"))

    # Total Hilbert space dimension
    dim = 2^n_qubits
    circuit = Matrix{ComplexF64}(I, dim, dim)

    for d in 1:depth
        if isodd(d)
            # Odd layers: pairs (1,2), (3,4), ...
            layer = Matrix{ComplexF64}(I, 4, 4)
            for q in 1:2:n_qubits
                if q < n_qubits
                    U2 = cue_matrix(4)
                    if q == 1
                        layer = U2
                    else
                        layer = kron(layer, U2)
                    end
                else
                    # Last unpaired qubit: identity
                    layer = kron(layer, Matrix{ComplexF64}(I, 2, 2))
                end
            end
        else
            # Even layers: pairs (2,3), (4,5), ...
            layer = Matrix{ComplexF64}(I, 2, 2)
            for q in 2:2:n_qubits
                if q < n_qubits
                    U2 = cue_matrix(4)
                    layer = kron(layer, U2)
                elseif q == n_qubits
                    # Last unpaired qubit: identity
                    layer = kron(layer, Matrix{ComplexF64}(I, 2, 2))
                end
            end
        end

        if d == 1
            circuit = layer
        else
            circuit = layer * circuit
        end
    end

    return circuit
end

end # module RandomUnitaries

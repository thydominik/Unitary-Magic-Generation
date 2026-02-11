"""
    circuits

Circuit construction and composition utilities.

This module focuses on building unitary matrices that represent circuit layouts.
Unitary sampling (e.g., Haar/CUE) is expected to come from `random_unitaries`, but this
module is written so that callers can inject their own two-qubit gate sampler.

All identifiers and docstrings are ASCII-only.
"""
module circuits

using LinearAlgebra
using Random

export brickwall_layer_matrix,
    brickwall_unitary_matrix

"""Return 2x2 identity as ComplexF64."""
function _i2()::Matrix{ComplexF64}
    return Matrix{ComplexF64}(I, 2, 2)
end

"""Fallback Haar sampler (used only if a sampler is not provided)."""
function _haar_unitary(dim::Integer; rng::AbstractRNG=Random.default_rng())::Matrix{ComplexF64}
    dim > 0 || throw(ArgumentError("dim must be positive"))

    z_re = randn(rng, dim, dim)
    z_im = randn(rng, dim, dim)
    z = ComplexF64.(z_re, z_im)

    f = qr(z)
    q = Matrix(f.Q)
    r = f.R

    d = diag(r)
    ph = similar(d)
    for i in eachindex(d)
        a = abs(d[i])
        ph[i] = a == 0 ? (1.0 + 0im) : (d[i] / a)
    end

    return q * Diagonal(ph)
end

"""Default 4x4 two-qubit gate sampler."""
function _default_two_qubit_gate(rng::AbstractRNG)::Matrix{ComplexF64}
    # Prefer random_unitaries.cue_matrix if the module is already loaded.
    if isdefined(Main, :random_unitaries)
        ru = getfield(Main, :random_unitaries)
        if ru isa Module && isdefined(ru, :cue_matrix)
            return getfield(ru, :cue_matrix)(4; rng=rng)
        end
    end

    return _haar_unitary(4; rng=rng)
end

"""
    brickwall_layer_matrix(n_qubits, layer_index; two_qubit_gate_sampler=_default_two_qubit_gate, rng=Random.default_rng())

Construct a single brick-wall layer as a unitary matrix on `n_qubits` qubits.

Layer pattern
- Odd layers (layer_index = 1, 3, ...): act on pairs (1,2), (3,4), ...
- Even layers (layer_index = 2, 4, ...): act on pairs (2,3), (4,5), ...

Open boundary conditions
- If the last qubit is unpaired in a given layer, an identity is applied on it.

Parameters
- `n_qubits`: Number of qubits.
- `layer_index`: Positive integer selecting odd/even pattern.
- `two_qubit_gate_sampler`: Function that takes an RNG and returns a 4x4 unitary.
- `rng`: Random-number generator.

Returns
- A unitary matrix of size `2^n_qubits x 2^n_qubits`.
"""
function brickwall_layer_matrix(
    n_qubits::Integer,
    layer_index::Integer;
    two_qubit_gate_sampler::Function=_default_two_qubit_gate,
    rng::AbstractRNG=Random.default_rng(),
)::Matrix{ComplexF64}
    n_qubits > 0 || throw(ArgumentError("n_qubits must be positive"))
    layer_index > 0 || throw(ArgumentError("layer_index must be positive"))

    i2 = _i2()

    blocks = Vector{Matrix{ComplexF64}}()

    if isodd(layer_index)
        q = 1
    else
        # Qubit 1 is unpaired in even layers.
        push!(blocks, i2)
        q = 2
    end

    while q <= n_qubits - 1
        push!(blocks, two_qubit_gate_sampler(rng))
        q += 2
    end

    if q == n_qubits
        push!(blocks, i2)
    end

    layer = blocks[1]
    for b in blocks[2:end]
        layer = kron(layer, b)
    end

    return layer
end

"""
    brickwall_unitary_matrix(n_qubits; depth=5*n_qubits, two_qubit_gate_sampler=_default_two_qubit_gate, rng=Random.default_rng())

Generate a brick-wall random unitary circuit as a full unitary matrix.

Parameters
- `n_qubits`: Number of qubits.
- `depth`: Number of layers.
- `two_qubit_gate_sampler`: Function that takes an RNG and returns a 4x4 unitary.
- `rng`: Random-number generator.

Returns
- A unitary matrix of size `2^n_qubits x 2^n_qubits`.

Notes
- This constructs the full matrix representation and scales exponentially in `n_qubits`.
"""
function brickwall_unitary_matrix(
    n_qubits::Integer;
    depth::Integer=5 * n_qubits,
    two_qubit_gate_sampler::Function=_default_two_qubit_gate,
    rng::AbstractRNG=Random.default_rng(),
)::Matrix{ComplexF64}
    n_qubits > 0 || throw(ArgumentError("n_qubits must be positive"))
    depth > 0 || throw(ArgumentError("depth must be positive"))

    dim = 2^n_qubits
    circuit = Matrix{ComplexF64}(I, dim, dim)

    for d in 1:depth
        layer = brickwall_layer_matrix(
            n_qubits,
            d;
            two_qubit_gate_sampler=two_qubit_gate_sampler,
            rng=rng,
        )
        circuit = layer * circuit
    end

    return circuit
end

end # module circuits

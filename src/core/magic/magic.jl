"""
    magic

Core routines to quantify stabilizer ("magic") properties of quantum states.

This module implements utilities to:
- Generate Pauli strings on `n` qubits.
- Build tensor-product Pauli operators.
- Compute magic measures for pure and mixed states based on stabiliser Rényi entropy.

The formulas follow the stabiliser Rényi entropy framework introduced in:
- https://doi.org/10.1103/PhysRevLett.128.050402

Notes
- The code is written for clarity and reproducibility and is suitable for research workflows.
- For performance-sensitive workloads, consider caching Pauli operator lists and avoiding repeated allocations.
"""
module magic

using LinearAlgebra
using SparseArrays

export generate_all_pauli_strings,
    pauli_matrix,
    pauli_operator_list,
    measure_magic_pure,
    measure_magic_mixed

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

"""Return `true` if `n` is a positive power of two."""
function _is_power_of_two(n::Integer)::Bool
    n > 0 && (n & (n - 1) == 0)
end

"""Infer number of qubits from a state-vector length; throws on invalid lengths."""
function _n_qubits_from_state_length(len::Integer)::Int
    _is_power_of_two(len) || throw(ArgumentError("state length must be a power of two, got $len"))
    return Int(round(log2(len)))
end

"""Infer number of qubits from a density-matrix size; throws on invalid sizes."""
function _n_qubits_from_density_matrix(ρ::AbstractMatrix)::Int
    size(ρ, 1) == size(ρ, 2) || throw(ArgumentError("density matrix must be square, got size $(size(ρ))"))
    _is_power_of_two(size(ρ, 1)) || throw(ArgumentError("density matrix dimension must be a power of two, got $(size(ρ, 1))"))
    return Int(round(log2(size(ρ, 1))))
end

# Precomputed single-qubit Pauli matrices (sparse), returned by reference.
const _PAULI_I = sparse(ComplexF64[1 0; 0 1])
const _PAULI_X = sparse(ComplexF64[0 1; 1 0])
const _PAULI_Y = sparse(ComplexF64[0 -im; im 0])
const _PAULI_Z = sparse(ComplexF64[1 0; 0 -1])

# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------

"""
    generate_all_pauli_strings(n_qubits::Integer) -> Vector{String}

Generate all Pauli strings on `n_qubits` qubits using the alphabet `{"I","X","Y","Z"}`.

Parameters
- `n_qubits`: Number of qubits (must be positive).

Returns
- A vector of strings of length `n_qubits`. The order matches `Iterators.product` over the
  per-qubit alphabet `["I","X","Y","Z"]`.
"""
function generate_all_pauli_strings(n_qubits::Integer)::Vector{String}
    n_qubits > 0 || throw(ArgumentError("n_qubits must be positive, got $n_qubits"))
    pauli_letters = ("I", "X", "Y", "Z")
    return [join(s) for s in Iterators.product(fill(pauli_letters, n_qubits)...)]
end

"""
    pauli_matrix(which) -> SparseMatrixCSC{ComplexF64,Int}

Return the single-qubit Pauli matrix specified by `which`.

Accepted labels
- `'I'`, `"I"`, or `0`
- `'X'`, `"X"`, or `1`
- `'Y'`, `"Y"`, or `2`
- `'Z'`, `"Z"`, or `3`

Returns
- A 2×2 sparse matrix of element type `ComplexF64`.
"""
function pauli_matrix(which)::SparseMatrixCSC{ComplexF64,Int}
    if which isa AbstractString
        length(which) == 1 || throw(ArgumentError("Pauli label string must have length 1, got '$which'"))
        which = only(which)
    end

    if which == 'I' || which == 0
        return _PAULI_I
    elseif which == 'X' || which == 1
        return _PAULI_X
    elseif which == 'Y' || which == 2
        return _PAULI_Y
    elseif which == 'Z' || which == 3
        return _PAULI_Z
    end

    throw(ArgumentError("Unsupported Pauli label: $which"))
end

"""
    pauli_operator_list(strings, n_qubits) -> Vector{SparseMatrixCSC{ComplexF64,Int}}

Construct tensor-product Pauli operators from Pauli strings.

Parameters
- `strings`: A vector of Pauli strings; each string must have length `n_qubits` and contain only
  `I`, `X`, `Y`, `Z`.
- `n_qubits`: Number of qubits (must be positive).

Returns
- A vector of sparse matrices of size `2^n_qubits × 2^n_qubits`.

Notes
- This function allocates `length(strings)` matrices and can be expensive for large `n_qubits`.
  Cache the result when reusing the same Pauli basis.
"""
function pauli_operator_list(strings::AbstractVector{<:AbstractString}, n_qubits::Integer)
    n_qubits > 0 || throw(ArgumentError("n_qubits must be positive, got $n_qubits"))

    ops = Vector{SparseMatrixCSC{ComplexF64,Int}}(undef, length(strings))
    for (i, s) in enumerate(strings)
        length(s) == n_qubits || throw(ArgumentError(
            "Pauli string length $(length(s)) does not match n_qubits=$n_qubits",
        ))

        op = pauli_matrix(s[1])
        for q in 2:n_qubits
            op = kron(op, pauli_matrix(s[q]))
        end
        ops[i] = sparse(op)
    end

    return ops
end

"""
    measure_magic_pure(state, pauli_ops, α=2) -> (magic, ξ)

Compute magic of a pure state using the stabiliser Rényi entropy construction.

Parameters
- `state`: State vector of length `2^n`.
- `pauli_ops`: Vector of Pauli operators on `n` qubits (e.g., from `pauli_operator_list`).
- `α`: Rényi parameter. Can be a single real number or a vector of real numbers.

Returns
- `magic`: A scalar if `α` is scalar, otherwise a vector.
- `ξ`: A vector with components
  \(ξ_P = \frac{1}{2^n} (\langle ψ | P | ψ \rangle)^2\) (real-valued as returned by the implementation).

Notes
- This function assumes the caller provides a properly normalized `state`.
- For `α == 1`, the implementation uses the special-case expression used in the legacy code.
"""
function measure_magic_pure(
    state::AbstractVector{<:Complex},
    pauli_ops::AbstractVector{<:AbstractMatrix},
    α::Union{Real, AbstractVector{<:Real}}=2,
)
    n_qubits = _n_qubits_from_state_length(length(state))
    hilbert_dim = 2^n_qubits

    ξ = Vector{Float64}(undef, length(pauli_ops))

    # Ensure we operate on a concrete vector to avoid surprises with views.
    ψ = ComplexF64.(state)
    ψ† = adjoint(ψ)

    for (i, op) in enumerate(pauli_ops)
        # Legacy convention: ξ ∝ (⟨ψ|P|ψ⟩)^2 (not abs2).
        v = (ψ† * (op * ψ))[1]
        ξ[i] = real((1 / hilbert_dim) * (v^2))
    end

    if α isa Real
        if α == 1
            magic = 1 - hilbert_dim * sum(ξ .^ 2)
        else
            a = float(α)
            magic = (1 - a)^(-1) * log2(sum(ξ .^ a)) - log2(hilbert_dim)
        end
    else
        magic = [(1 - float(a))^(-1) * log2(sum(ξ .^ float(a))) - log2(hilbert_dim) for a in α]
    end

    return magic, ξ
end

"""
    measure_magic_mixed(ρ, pauli_ops; α=2) -> Float64

Compute magic of a mixed state `ρ` using stabiliser Rényi entropy.

Parameters
- `ρ`: Density matrix of size `2^n × 2^n`.
- `pauli_ops`: Vector of Pauli operators on `n` qubits.
- `α`: Rényi parameter (keyword; default `2`).

Returns
- Magic value as `Float64`.

Notes
- This implementation follows the legacy code structure and includes the purity correction term
  `log2(sum(eigvals(ρ).^2))`.
- For large systems, eigenvalue computations may dominate runtime.
"""
function measure_magic_mixed(
    ρ::AbstractMatrix{<:Complex},
    pauli_ops::AbstractVector{<:AbstractMatrix};
    α::Real=2,
)::Float64
    n_qubits = _n_qubits_from_density_matrix(ρ)
    hilbert_dim = 2^n_qubits

    ξ = Vector{Float64}(undef, length(pauli_ops))
    for (i, op) in enumerate(pauli_ops)
        # trace(op * ρ) implemented via sum(diag(...)) for broad compatibility.
        t = sum(diag(op * ρ))
        ξ[i] = real((1 / hilbert_dim) * (t^2))
    end

    magic = (1 - α)^(-1) * log2(sum(ξ .^ α)) - log2(hilbert_dim) + log2(sum(eigvals(ρ) .^ 2))
    return magic
end

end # module magic

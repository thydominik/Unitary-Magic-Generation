"""
    magic

Core routines to quantify stabilizer ("magic") properties of quantum states.

This module implements utilities to:
- Generate Pauli strings on `n` qubits.
- Build tensor-product Pauli operators.
- Compute magic measures for pure and mixed states based on stabiliser Renyi entropy.

The formulas follow the stabiliser Renyi entropy framework introduced in:
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
function _n_qubits_from_density_matrix(rho::AbstractMatrix)::Int
    size(rho, 1) == size(rho, 2) || throw(ArgumentError("density matrix must be square, got size $(size(rho))"))
    _is_power_of_two(size(rho, 1)) || throw(ArgumentError("density matrix dimension must be a power of two, got $(size(rho, 1))"))
    return Int(round(log2(size(rho, 1))))
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

Generate all Pauli strings on `n_qubits` qubits using the alphabet `{I,X,Y,Z}`.

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
- A 2x2 sparse matrix of element type `ComplexF64`.
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
- A vector of sparse matrices of size `2^n_qubits x 2^n_qubits`.

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
    measure_magic_pure(state, pauli_ops, alpha=2) -> (magic, xi)

Compute magic of a pure state using the stabiliser Renyi entropy construction.

Parameters
- `state`: State vector of length `2^n`.
- `pauli_ops`: Vector of Pauli operators on `n` qubits (e.g., from `pauli_operator_list`).
- `alpha`: Renyi parameter. Can be a single real number or a vector of real numbers.

Returns
- `magic`: A scalar if `alpha` is scalar, otherwise a vector.
- `xi`: A vector with components
  `xi_P = (1 / 2^n) * (<psi|P|psi>)^2` (real-valued as returned by the implementation).

Notes
- This function assumes the caller provides a properly normalized `state`.
- For `alpha == 1`, the implementation uses the special-case expression used in the legacy code.
"""
function measure_magic_pure(
    state::AbstractVector{<:Complex},
    pauli_ops::AbstractVector{<:AbstractMatrix},
    alpha::Union{Real, AbstractVector{<:Real}}=2,
)
    n_qubits = _n_qubits_from_state_length(length(state))
    hilbert_dim = 2^n_qubits

    xi = Vector{Float64}(undef, length(pauli_ops))

    # Ensure we operate on a concrete vector to avoid surprises with views.
    psi = ComplexF64.(state)
    psi_dag = adjoint(psi)

    for (i, op) in enumerate(pauli_ops)
        # Legacy convention: xi is proportional to (<psi|P|psi>)^2 (not abs2).
        v = (psi_dag * (op * psi))[1]
        xi[i] = real((1 / hilbert_dim) * (v^2))
    end

    if alpha isa Real
        if alpha == 1
            magic = 1 - hilbert_dim * sum(xi .^ 2)
        else
            a = float(alpha)
            magic = (1 - a)^(-1) * log2(sum(xi .^ a)) - log2(hilbert_dim)
        end
    else
        magic = [(1 - float(a))^(-1) * log2(sum(xi .^ float(a))) - log2(hilbert_dim) for a in alpha]
    end

    return magic, xi
end

"""
    measure_magic_mixed(rho, pauli_ops; alpha=2) -> Float64

Compute magic of a mixed state `rho` using stabiliser Renyi entropy.

Parameters
- `rho`: Density matrix of size `2^n x 2^n`.
- `pauli_ops`: Vector of Pauli operators on `n` qubits.
- `alpha`: Renyi parameter (keyword; default `2`).

Returns
- Magic value as `Float64`.

Notes
- This implementation follows the legacy code structure and includes the purity correction term
  `log2(sum(eigvals(rho).^2))`.
- For large systems, eigenvalue computations may dominate runtime.
"""
function measure_magic_mixed(
    rho::AbstractMatrix{<:Complex},
    pauli_ops::AbstractVector{<:AbstractMatrix};
    alpha::Real=2,
)::Float64
    n_qubits = _n_qubits_from_density_matrix(rho)
    hilbert_dim = 2^n_qubits

    xi = Vector{Float64}(undef, length(pauli_ops))
    for (i, op) in enumerate(pauli_ops)
        # trace(op * rho) implemented via sum(diag(...)) for broad compatibility.
        t = sum(diag(op * rho))
        xi[i] = real((1 / hilbert_dim) * (t^2))
    end

    magic = (1 - alpha)^(-1) * log2(sum(xi .^ alpha)) - log2(hilbert_dim) + log2(sum(eigvals(rho) .^ 2))
    return magic
end

end # module magic

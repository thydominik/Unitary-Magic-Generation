"""
    entanglement

Core routines for bipartite entanglement and entropy quantification of pure quantum states.

This module provides:
- Reduced density matrices via partial trace (pure-state optimized implementation).
- von Neumann entropy.
- Renyi entropies.
- Convenience functions to evaluate entanglement/entropy across standard bipartitions.

Notes
- All identifiers and docstrings are ASCII-only by design.
- Qubit indexing is 1-based. A state vector of length 2^n is interpreted in the computational
  basis with qubit 1 as the least significant bit, consistent with Julia's column-major ordering
  and the legacy implementation in this repository.
"""
module entanglement

using LinearAlgebra

export reduced_density_matrix,
    von_neumann_entropy,
    renyi_entropy,
    calculate_entanglement,
    calculate_renyi_entropy

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

"""Return true if `n` is a positive power of two."""
function _is_power_of_two(n::Integer)::Bool
    n > 0 && (n & (n - 1) == 0)
end

"""Infer number of qubits from a state-vector length; throws on invalid lengths."""
function _n_qubits_from_state_length(len::Integer)::Int
    _is_power_of_two(len) || throw(ArgumentError("state length must be a power of two, got $len"))
    return Int(round(log2(len)))
end

"""Validate qubit indices and return a sorted, unique copy."""
function _validate_qubit_indices(indices::AbstractVector{<:Integer}, n_qubits::Integer)::Vector{Int}
    all(1 .<= indices .<= n_qubits) || throw(ArgumentError("qubit indices must be between 1 and $n_qubits"))
    unique_indices = unique(Int.(indices))
    length(unique_indices) == length(indices) || throw(ArgumentError("qubit indices must be unique"))
    sort!(unique_indices)
    return unique_indices
end

"""Return eigenvalues filtered to remove tiny negative/complex numerical noise."""
function _filtered_eigenvalues(rho::AbstractMatrix{<:Complex}; tol::Real=1e-14)
    vals = eigvals(rho)
    # Keep values with positive real part; take real part for stability.
    reals = Float64[]
    for v in vals
        r = real(v)
        if r > tol
            push!(reals, r)
        end
    end
    return reals
end

# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------

"""
    reduced_density_matrix(state_vector, subsystem_qubits) -> (rho_reduced, keep_qubits)

Compute the reduced density matrix by tracing out `subsystem_qubits` from a pure state.

Parameters
- `state_vector`: Pure state vector of length 2^n.
- `subsystem_qubits`: Indices (1-based) of qubits to trace out.

Returns
- `rho_reduced`: Reduced density matrix on the remaining qubits.
- `keep_qubits`: The qubits that were kept (complement of `subsystem_qubits`).

Implementation
This function avoids constructing the full 2^n x 2^n density matrix. It reshapes the state
into a tensor, permutes axes to group kept/traced qubits, reshapes into a matrix Psi, and
computes `rho_reduced = Psi * Psi'`.
"""
function reduced_density_matrix(
    state_vector::AbstractVector{<:Complex},
    subsystem_qubits::AbstractVector{<:Integer},
)
    n_qubits = _n_qubits_from_state_length(length(state_vector))

    traced = _validate_qubit_indices(subsystem_qubits, n_qubits)
    keep = setdiff(collect(1:n_qubits), traced)

    # Special case: trace out nothing.
    if isempty(traced)
        psi = ComplexF64.(state_vector)
        rho = psi * adjoint(psi)
        return rho, keep
    end

    # Special case: trace out everything.
    if isempty(keep)
        return reshape(ComplexF64[1.0 + 0im], 1, 1), keep
    end

    dims = ntuple(_ -> 2, n_qubits)
    psi_tensor = reshape(ComplexF64.(state_vector), dims)

    perm = vcat(keep, traced)
    psi_perm = permutedims(psi_tensor, perm)

    dim_keep = 2^length(keep)
    dim_traced = 2^length(traced)
    psi_mat = reshape(psi_perm, dim_keep, dim_traced)

    rho_reduced = psi_mat * adjoint(psi_mat)
    return rho_reduced, keep
end

"""
    von_neumann_entropy(rho) -> Float64

Compute von Neumann entropy S(rho) = -Tr(rho log2(rho)).

Parameters
- `rho`: Density matrix (square). This function uses eigenvalues and ignores tiny numerical
  eigenvalues below a fixed tolerance.

Returns
- Entropy in bits.
"""
function von_neumann_entropy(rho::AbstractMatrix{<:Complex})::Float64
    vals = _filtered_eigenvalues(rho)
    if isempty(vals)
        return 0.0
    end
    return -sum(p -> p * log2(p), vals)
end

"""
    renyi_entropy(rho, alpha) -> Float64

Compute Renyi entropy of order `alpha` for a density matrix.

Definition
- For alpha == 1: returns von Neumann entropy.
- For alpha == 0: returns log2(rank(rho)) based on eigenvalues above tolerance.
- Otherwise: S_alpha(rho) = 1/(1-alpha) * log2(sum_i lambda_i^alpha).

Parameters
- `rho`: Density matrix.
- `alpha`: Renyi order (must be non-negative).
"""
function renyi_entropy(rho::AbstractMatrix{<:Complex}, alpha::Real)::Float64
    alpha >= 0 || throw(ArgumentError("Renyi index alpha must be non-negative"))

    # alpha == 1 (limit) -> von Neumann entropy
    if abs(alpha - 1.0) < 1e-10
        return von_neumann_entropy(rho)
    end

    vals = _filtered_eigenvalues(rho)
    isempty(vals) && return 0.0

    # alpha == 0 -> log2(rank)
    if abs(alpha) < 1e-12
        return log2(length(vals))
    end

    s = sum(v -> v^alpha, vals)
    return (1 / (1 - alpha)) * log2(s)
end

"""
    calculate_entanglement(state_vector; subsystems="All") -> (entropies, kept_qubits)

Compute von Neumann entanglement entropies across standard bipartitions of a pure state.

Parameters
- `state_vector`: Pure state vector of length 2^n.
- `subsystems`: "All" to evaluate cuts tracing out qubits 1:k for k=1..n-1, or "Middle" for the
  single cut tracing out qubits 1:floor(n/2).

Returns
- `entropies`: Vector of entropies for each cut.
- `kept_qubits`: Vector of kept-qubit index lists corresponding to each entropy.
"""
function calculate_entanglement(
    state_vector::AbstractVector{<:Complex};
    subsystems::String="All",
)
    n_qubits = _n_qubits_from_state_length(length(state_vector))
    n_qubits >= 2 || throw(ArgumentError("need at least 2 qubits to compute bipartite entanglement"))
    subsystems in ("All", "Middle") || throw(ArgumentError("subsystems must be 'All' or 'Middle'"))

    entropies = Float64[]
    kept = Vector{Vector{Int}}()

    if subsystems == "All"
        for k in 1:(n_qubits - 1)
            rho_red, keep_qubits = reduced_density_matrix(state_vector, 1:k)
            push!(entropies, von_neumann_entropy(rho_red))
            push!(kept, keep_qubits)
        end
    else
        k = floor(Int, n_qubits / 2)
        rho_red, keep_qubits = reduced_density_matrix(state_vector, 1:k)
        push!(entropies, von_neumann_entropy(rho_red))
        push!(kept, keep_qubits)
    end

    return entropies, kept
end

"""
    calculate_renyi_entropy(state_vector, alpha; subsystems="All") -> Vector{Float64}

Compute Renyi entropies of order `alpha` across standard bipartitions of a pure state.

Parameters
- `state_vector`: Pure state vector of length 2^n.
- `alpha`: Renyi order.
- `subsystems`: "All" or "Middle" (same meaning as in `calculate_entanglement`).

Returns
- Vector of Renyi entropies for each cut.
"""
function calculate_renyi_entropy(
    state_vector::AbstractVector{<:Complex},
    alpha::Real;
    subsystems::String="All",
)::Vector{Float64}
    n_qubits = _n_qubits_from_state_length(length(state_vector))
    n_qubits >= 2 || throw(ArgumentError("need at least 2 qubits to compute bipartite entropies"))
    subsystems in ("All", "Middle") || throw(ArgumentError("subsystems must be 'All' or 'Middle'"))
    alpha >= 0 || throw(ArgumentError("Renyi index alpha must be non-negative"))

    entropies = Float64[]

    if subsystems == "All"
        for k in 1:(n_qubits - 1)
            rho_red, _ = reduced_density_matrix(state_vector, 1:k)
            push!(entropies, renyi_entropy(rho_red, alpha))
        end
    else
        k = floor(Int, n_qubits / 2)
        rho_red, _ = reduced_density_matrix(state_vector, 1:k)
        push!(entropies, renyi_entropy(rho_red, alpha))
    end

    return entropies
end

end # module entanglement

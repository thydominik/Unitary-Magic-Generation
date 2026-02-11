"""
    entanglement

Core routines for bipartite entanglement and entropy quantification of pure quantum states.

This module provides:
- Reduced density matrices via partial trace (pure-state optimized implementation).
- von Neumann entropy.
- Renyi entropies.
- Entanglement negativity and logarithmic negativity (pure-state optimized via Schmidt coefficients).
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
    entanglement_negativity,
    logarithmic_negativity,
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

"""
    _schmidt_coefficients(state_vector, subsystem_qubits) -> Vector{Float64}

Compute Schmidt coefficients for the bipartition defined by `subsystem_qubits` and its complement.

This uses a reshape + permute + SVD approach and is appropriate for pure states.
"""
function _schmidt_coefficients(
    state_vector::AbstractVector{<:Complex},
    subsystem_qubits::AbstractVector{<:Integer},
)::Vector{Float64}
    n_qubits = _n_qubits_from_state_length(length(state_vector))

    a = _validate_qubit_indices(subsystem_qubits, n_qubits)
    b = setdiff(collect(1:n_qubits), a)

    # If one side is empty, there is no entanglement.
    if isempty(a) || isempty(b)
        return Float64[1.0]
    end

    psi = ComplexF64.(state_vector)
    nrm = norm(psi)
    nrm == 0 && throw(ArgumentError("state vector must be nonzero"))
    psi ./= nrm

    dims = ntuple(_ -> 2, n_qubits)
    psi_tensor = reshape(psi, dims)

    # Group qubits as [b, a]. Negativity is symmetric, so either order is fine.
    perm = vcat(b, a)
    psi_perm = permutedims(psi_tensor, perm)

    dim_b = 2^length(b)
    dim_a = 2^length(a)
    psi_mat = reshape(psi_perm, dim_b, dim_a)

    s = svdvals(psi_mat)
    return Float64.(s)
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
    entanglement_negativity(state_vector, subsystem_qubits) -> Float64

Compute the entanglement negativity for a pure state across a bipartition.

This function returns the (non-logarithmic) negativity

    N = (||rho_PT||_1 - 1) / 2

where rho_PT is the partial transpose of the bipartite density matrix. For pure states,
this can be computed from the Schmidt coefficients s_i as

    N = ((sum_i s_i)^2 - 1) / 2.

Parameters
- `state_vector`: Pure state vector of length 2^n.
- `subsystem_qubits`: Qubit indices that define one side of the bipartition. The result is
  symmetric under swapping the subsystem and its complement.

Returns
- Negativity as Float64.
"""
function entanglement_negativity(
    state_vector::AbstractVector{<:Complex},
    subsystem_qubits::AbstractVector{<:Integer},
)::Float64
    s = _schmidt_coefficients(state_vector, subsystem_qubits)
    sum_s = sum(s)

    # Numerical safety: sum_s should be >= 1 for normalized pure states.
    if sum_s < 1.0
        sum_s = 1.0
    end

    return 0.5 * (sum_s^2 - 1.0)
end

"""
    logarithmic_negativity(state_vector, subsystem_qubits) -> Float64

Compute the logarithmic entanglement negativity for a pure state across a bipartition.

Definition
- E_N = log2(||rho_PT||_1)

For pure states, if s_i are Schmidt coefficients then
- ||rho_PT||_1 = (sum_i s_i)^2
- E_N = 2 * log2(sum_i s_i)

Parameters
- `state_vector`: Pure state vector of length 2^n.
- `subsystem_qubits`: Qubit indices that define one side of the bipartition.

Returns
- Logarithmic negativity in bits.
"""
function logarithmic_negativity(
    state_vector::AbstractVector{<:Complex},
    subsystem_qubits::AbstractVector{<:Integer},
)::Float64
    s = _schmidt_coefficients(state_vector, subsystem_qubits)
    sum_s = sum(s)

    if sum_s < 1.0
        sum_s = 1.0
    end

    return 2.0 * log2(sum_s)
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

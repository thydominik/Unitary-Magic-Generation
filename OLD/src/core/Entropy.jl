"""
    Entropy

Rényi entropy quantification for quantum states.

This module provides tools for computing Rényi entropies of various orders across
bipartite cuts of pure quantum states. The von Neumann entropy emerges as the
\\(\\alpha \\to 1\\) limit.

# Main Functions

- [`renyi_entropy`](@ref): Compute Rényi entropy of order \\(\\alpha\\) for a density matrix
- [`calculate_renyi_entropy`](@ref): Compute Rényi entropy across bipartitions

# References

- Rényi, A., "On measures of entropy and information", Proc. 4th Berkeley Symp. Math. Stat. Prob. (1961)
- Müller-Lennert et al., "On quantum Rényi entropies: A new generalization and some properties", J. Math. Phys. 54, 122203 (2013)
"""
module Entropy

using LinearAlgebra

export renyi_entropy,
       calculate_renyi_entropy

# Import reduced_density_matrix from Entanglement module
using ..Entanglement: reduced_density_matrix, von_neumann_entropy

"""
    renyi_entropy(density_matrix::Matrix{<:Complex}, α::Real) -> Float64

Compute the Rényi entropy of order \\(\\alpha\\) for a density matrix.

The Rényi entropy is defined as:
\[
S_\\alpha(\\rho) = \\frac{1}{1 - \\alpha} \\log_2 \\left( \\sum_i \\lambda_i^\\alpha \\right)
\]
where \\(\\lambda_i\\) are the eigenvalues of \\(\\rho\\).

For \\(\\alpha = 1\\), this reduces to the von Neumann entropy:
\[
S_1(\\rho) = -\\sum_i \\lambda_i \\log_2 \\lambda_i
\]

# Arguments

- `density_matrix::Matrix{<:Complex}`: Density matrix (Hermitian, positive semi-definite)
- `α::Real`: Rényi index (order). Common values:
  - `α = 0`: Max-entropy (log of rank)
  - `α = 1`: von Neumann entropy (limiting case)
  - `α = 2`: Collision entropy
  - `α → ∞`: Min-entropy

# Returns

- `entropy::Float64`: Rényi entropy in bits (base-2 logarithm)

# Example

```julia
# Maximally mixed state on 1 qubit
ρ = [0.5 0.0; 0.0 0.5]
S2 = renyi_entropy(ρ, 2.0)  # Collision entropy, should return 1.0
```
"""
function renyi_entropy(density_matrix::Matrix{<:Complex}, α::Real)
    if α < 0
        throw(ArgumentError("Rényi index α must be non-negative"))
    end
    
    eigenvalues = eigvals(density_matrix)
    
    # Filter out negligible eigenvalues
    λ_filtered = filter(x -> real(x) > 1e-14, eigenvalues)
    
    if isempty(λ_filtered)
        return 0.0
    end
    
    # Special case: α = 1 corresponds to von Neumann entropy
    if abs(α - 1.0) < 1e-10
        return von_neumann_entropy(density_matrix)
    end
    
    # General Rényi entropy formula
    entropy = (1 / (1 - α)) * log2(sum(real(λ)^α for λ in λ_filtered))
    
    return entropy
end

"""
    calculate_renyi_entropy(state_vector::AbstractVector{<:Complex}, α::Real; subsystems::String = "All") -> Vector{Float64}

Calculate Rényi entropy of order \\(\\alpha\\) across bipartite cuts.

# Arguments

- `state_vector::AbstractVector{<:Complex}`: Pure state vector of length `2^n`
- `α::Real`: Rényi index (order)
- `subsystems::String = "All"`: Which bipartitions to compute
  - `"All"`: All cuts from 1:|n-1| qubits to |n-1|:1 qubits
  - `"Middle"`: Only the middle cut (⌊n/2⌋ qubits traced out)

# Returns

- `entropies::Vector{Float64}`: Rényi entropies for each bipartition

# Example

```julia
using LinearAlgebra

# 3-qubit W state |001⟩ + |010⟩ + |100⟩
ψ = zeros(ComplexF64, 8)
ψ[2] = ψ[3] = ψ[5] = 1/sqrt(3)

S2_all = calculate_renyi_entropy(ψ, 2.0; subsystems="All")
# Collision entropy across all bipartitions
```
"""
function calculate_renyi_entropy(
    state_vector::AbstractVector{<:Complex},
    α::Real;
    subsystems::String = "All"
)
    n_qubits = Int(log2(length(state_vector)))
    
    if 2^n_qubits != length(state_vector)
        throw(ArgumentError("State vector length must be a power of 2"))
    end
    if n_qubits < 2
        throw(ArgumentError("Need at least 2 qubits to compute entropy across bipartitions"))
    end
    if !(subsystems in ["All", "Middle"])
        throw(ArgumentError("subsystems must be \"All\" or \"Middle\""))
    end
    if α < 0
        throw(ArgumentError("Rényi index α must be non-negative"))
    end
    
    entropies = Float64[]
    
    if subsystems == "All"
        for k in 1:(n_qubits - 1)
            subsystem = 1:k
            reduced_matrix, _ = reduced_density_matrix(state_vector, collect(subsystem))
            push!(entropies, renyi_entropy(reduced_matrix, α))
        end
    elseif subsystems == "Middle"
        subsystem = 1:floor(Int, n_qubits / 2)
        reduced_matrix, _ = reduced_density_matrix(state_vector, collect(subsystem))
        push!(entropies, renyi_entropy(reduced_matrix, α))
    end
    
    return entropies
end

end # module Entropy

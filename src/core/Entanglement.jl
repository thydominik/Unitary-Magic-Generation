"""
    Entanglement

Entanglement quantification for quantum states.

This module provides tools for computing reduced density matrices and entanglement entropies
of pure quantum states across bipartite cuts.

# Main Functions

- [`reduced_density_matrix`](@ref): Compute reduced density matrix by tracing out subsystems
- [`von_neumann_entropy`](@ref): Calculate von Neumann entropy of a density matrix
- [`calculate_entanglement`](@ref): Compute entanglement entropy across multiple bipartitions

# References

- Nielsen & Chuang, "Quantum Computation and Quantum Information" (2000)
- Horodecki et al., "Quantum entanglement", Rev. Mod. Phys. 81, 865 (2009)
"""
module Entanglement

using LinearAlgebra
using Combinatorics

export reduced_density_matrix,
       von_neumann_entropy,
       calculate_entanglement

"""
    reduced_density_matrix(state_vector::AbstractVector{<:Complex}, subsystem_qubits::Vector{Int}) -> (Matrix{ComplexF64}, Vector{Int})

Compute the reduced density matrix by tracing out specified qubits.

# Arguments

- `state_vector::AbstractVector{<:Complex}`: Pure state vector of length `2^n` for `n` qubits
- `subsystem_qubits::Vector{Int}`: Indices of qubits to trace out (1-indexed)

# Returns

- `reduced_matrix::Matrix{ComplexF64}`: Reduced density matrix after partial trace
- `keep_qubits::Vector{Int}`: Indices of qubits kept in the reduced system

# Notes

The partial trace is performed by summing over all computational basis states of the
traced-out subsystem. The remaining qubits form the reduced density matrix.

# Example

```julia
using LinearAlgebra

# Bell state |Φ+⟩ = (|00⟩ + |11⟩)/√2
ψ = [1, 0, 0, 1] / sqrt(2)

# Trace out qubit 2
ρ_reduced, kept = reduced_density_matrix(ψ, [2])
# Result: maximally mixed state on qubit 1
```
"""
function reduced_density_matrix(
    state_vector::AbstractVector{<:Complex},
    subsystem_qubits::Vector{Int}
)
    n_qubits = Int(log2(length(state_vector)))
    
    # Input validation
    if 2^n_qubits != length(state_vector)
        throw(ArgumentError("State vector length must be a power of 2"))
    end
    if !all(1 .<= subsystem_qubits .<= n_qubits)
        throw(ArgumentError("Subsystem qubit indices must be between 1 and $n_qubits"))
    end
    if length(unique(subsystem_qubits)) != length(subsystem_qubits)
        throw(ArgumentError("Subsystem qubit indices must be unique"))
    end
    
    # Construct full density matrix
    full_density_matrix = state_vector * state_vector'
    
    # Qubits to keep in the reduced system
    keep_qubits = setdiff(1:n_qubits, subsystem_qubits)
    
    # Size of the reduced density matrix
    reduced_size = 2^length(keep_qubits)
    reduced_matrix = zeros(ComplexF64, reduced_size, reduced_size)
    
    # Perform partial trace
    for (i, basis_i) in enumerate(Iterators.product(fill(0:1, length(keep_qubits))...))
        for (j, basis_j) in enumerate(Iterators.product(fill(0:1, length(keep_qubits))...))
            sum_val = zero(ComplexF64)
            
            # Sum over traced-out subsystem basis states
            for traced_basis in Iterators.product(fill(0:1, length(subsystem_qubits))...)
                full_i = [basis_i..., traced_basis...]
                full_j = [basis_j..., traced_basis...]
                
                perm_i = sortperm([keep_qubits..., subsystem_qubits...])
                perm_j = sortperm([keep_qubits..., subsystem_qubits...])
                
                idx_i = sum(full_i[perm_i] .* (2 .^ (0:n_qubits-1))) + 1
                idx_j = sum(full_j[perm_j] .* (2 .^ (0:n_qubits-1))) + 1
                
                sum_val += full_density_matrix[idx_i, idx_j]
            end
            
            reduced_matrix[i, j] = sum_val
        end
    end
    
    return reduced_matrix, keep_qubits
end

"""
    von_neumann_entropy(density_matrix::Matrix{<:Complex}) -> Float64

Calculate the von Neumann entropy of a density matrix.

The von Neumann entropy is defined as:
\[
S(\\rho) = -\\text{Tr}(\\rho \\log_2 \\rho) = -\\sum_i \\lambda_i \\log_2 \\lambda_i
\]
where \\(\\lambda_i\\) are the eigenvalues of \\(\\rho\\).

# Arguments

- `density_matrix::Matrix{<:Complex}`: Density matrix (Hermitian, positive semi-definite)

# Returns

- `entropy::Float64`: von Neumann entropy in bits (base-2 logarithm)

# Example

```julia
# Maximally mixed state on 1 qubit
ρ = [0.5 0.0; 0.0 0.5]
S = von_neumann_entropy(ρ)  # Should return 1.0
```
"""
function von_neumann_entropy(density_matrix::Matrix{<:Complex})
    eigenvalues = eigvals(density_matrix)
    return -sum(λ -> real(λ) > 1e-14 ? real(λ) * log2(real(λ)) : 0.0, eigenvalues)
end

"""
    calculate_entanglement(state_vector::AbstractVector{<:Complex}; subsystems::String = "All") -> (Vector{Float64}, Vector{Vector{Int}})

Calculate entanglement entropy across bipartite cuts of a quantum state.

# Arguments

- `state_vector::AbstractVector{<:Complex}`: Pure state vector of length `2^n`
- `subsystems::String = "All"`: Which bipartitions to compute
  - `"All"`: All cuts from 1:|n-1| qubits to |n-1|:1 qubits
  - `"Middle"`: Only the middle cut (⌊n/2⌋ qubits traced out)

# Returns

- `entropies::Vector{Float64}`: von Neumann entropies for each bipartition
- `subsystems::Vector{Vector{Int}}`: Kept qubit indices for each bipartition

# Example

```julia
using LinearAlgebra

# 3-qubit GHZ state |000⟩ + |111⟩
ψ = zeros(ComplexF64, 8)
ψ[1] = ψ[8] = 1/sqrt(2)

entropies, subsystems = calculate_entanglement(ψ; subsystems="All")
# Returns entropies for 1:2 and 2:1 cuts
```
"""
function calculate_entanglement(
    state_vector::AbstractVector{<:Complex};
    subsystems::String = "All"
)
    n_qubits = Int(log2(length(state_vector)))
    
    if 2^n_qubits != length(state_vector)
        throw(ArgumentError("State vector length must be a power of 2"))
    end
    if n_qubits < 2
        throw(ArgumentError("Need at least 2 qubits to compute entanglement"))
    end
    if !(subsystems in ["All", "Middle"])
        throw(ArgumentError("subsystems must be \"All\" or \"Middle\""))
    end
    
    entropies = Float64[]
    kept_subsystems = Vector{Int}[]
    
    if subsystems == "All"
        for k in 1:(n_qubits - 1)
            subsystem = 1:k
            reduced_matrix, keep_qubits = reduced_density_matrix(state_vector, collect(subsystem))
            push!(entropies, von_neumann_entropy(reduced_matrix))
            push!(kept_subsystems, keep_qubits)
        end
    elseif subsystems == "Middle"
        subsystem = 1:floor(Int, n_qubits / 2)
        reduced_matrix, keep_qubits = reduced_density_matrix(state_vector, collect(subsystem))
        push!(entropies, von_neumann_entropy(reduced_matrix))
        push!(kept_subsystems, keep_qubits)
    end
    
    return entropies, kept_subsystems
end

end # module Entanglement

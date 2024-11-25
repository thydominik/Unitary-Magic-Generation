module Measure_Entanglement

using LinearAlgebra
using Combinatorics

export reduced_density_matrix, von_neumann_entropy, calculate_entanglement

function reduced_density_matrix(state_vector::Vector{ComplexF64}, subsystem_qubits::Vector{Int})

    No_Qubits = Int(log2(length(state_vector)))

    # Costruct the full density matrix using the state input
    full_density_matrix = state_vector * state_vector'

    #The qubits that we keep as a subsystem:
    keep_qubits = setdiff(1:No_Qubits, subsystem_qubits)
    
    # Size of the reduced density matrix
    reduced_size = 2^length(keep_qubits)
    # Space allocation
    reduced_matrix = zeros(Complex{Float64}, reduced_size, reduced_size)
    
    for (i, basis_i) in enumerate(Iterators.product(fill(0:1, length(keep_qubits))...))
        for (j, basis_j) in enumerate(Iterators.product(fill(0:1, length(keep_qubits))...))
            sum_val = 0.0 + 0.0im
            for traced_basis in Iterators.product(fill(0:1, length(subsystem_qubits))...)
                full_i = [basis_i..., traced_basis...]
                full_j = [basis_j..., traced_basis...]
                perm_i = sortperm([keep_qubits..., subsystem_qubits...])
                perm_j = sortperm([keep_qubits..., subsystem_qubits...])
                idx_i = sum(full_i[perm_i] .* (2 .^ (0:No_Qubits-1))) + 1
                idx_j = sum(full_j[perm_j] .* (2 .^ (0:No_Qubits-1))) + 1
                sum_val += full_density_matrix[idx_i, idx_j]
            end
            reduced_matrix[i, j] = sum_val
        end
    end
    return reduced_matrix, keep_qubits
end

function von_neumann_entropy(density_matrix::Matrix{Complex{Float64}})
    eigenvalues = eigvals(density_matrix)
    return -sum(λ -> λ > 0 ? λ * log2(λ) : 0, eigenvalues)
end

function calculate_entanglement(state_vector::Vector{Complex{Float64}})
    No_Qubits   = Int(log2(length(state_vector)))
    Entropies   = Float64[]
    SubSystems  = Vector{Int}[]

    for k in 1:(No_Qubits - 1)
        subsystem = 1:k
        reduced_matrix, keep_qubits = reduced_density_matrix(state_vector, collect(subsystem))
        push!(Entropies, von_neumann_entropy(reduced_matrix))
        push!(SubSystems, keep_qubits)
    end
    
    return Entropies, SubSystems
end

end
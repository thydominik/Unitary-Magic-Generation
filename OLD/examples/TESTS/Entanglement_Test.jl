using LinearAlgebra

function reduced_density_matrix(state_vector::Vector{Complex{T}}, subsystem_size::Int) where T
    # Calculate the number of qubits
    N = Int(log2(length(state_vector)))
    
    # Create the full density matrix from the state vector
    full_density_matrix = state_vector * transpose(conj(state_vector))
    
    # Define the dimension of the reduced matrix
    reduced_size = 2^(N - subsystem_size)
    reduced_matrix = zeros(Complex{T}, reduced_size, reduced_size)
    
    # Tracing out the subsystem
    for i in 0:(reduced_size - 1)
        for j in 0:(reduced_size - 1)
            for k in 0:(2^subsystem_size - 1)
                index_i = (i << subsystem_size) | k
                index_j = (j << subsystem_size) | k
                reduced_matrix[i + 1, j + 1] += full_density_matrix[index_i + 1, index_j + 1]
            end
        end
    end
    
    return reduced_matrix
end

function von_neumann_entropy(reduced_matrix::Matrix{Complex{T}}) where T
    # Calculate eigenvalues of the reduced density matrix
    eigenvalues = eigen(reduced_matrix).values
    
    # Filter out zero eigenvalues and calculate entropy
    entropy = -sum(eigenvalue * log2(eigenvalue) for eigenvalue in eigenvalues if eigenvalue > 0)
    
    return entropy
end


N = 5
Ψ = rand(ComplexF64, 2^N)
Ψ = Ψ ./ norm(Ψ)
E = []

for Bip in 1:N-1
    ϱ = reduced_density_matrix(Ψ, Bip)
    push!(E, von_neumann_entropy(ϱ))
    
end

println(mean(E))









using LinearAlgebra
using Random
using Combinatorics
using Statistics
# Function to calculate reduced density matrix

function reduced_density_matrix(state_vector::Vector{Complex{Float64}}, subsystem_qubits::Vector{Int})
    N = Int(log2(length(state_vector)))
    full_density = state_vector * state_vector'
    keep_qubits = setdiff(1:N, subsystem_qubits)
    
    reduced_size = 2^length(keep_qubits)
    reduced_matrix = zeros(Complex{Float64}, reduced_size, reduced_size)
    
    for (i, basis_i) in enumerate(Iterators.product(fill(0:1, length(keep_qubits))...))
        for (j, basis_j) in enumerate(Iterators.product(fill(0:1, length(keep_qubits))...))
            sum_val = 0.0 + 0.0im
            for traced_basis in Iterators.product(fill(0:1, length(subsystem_qubits))...)
                full_i = [basis_i..., traced_basis...]
                full_j = [basis_j..., traced_basis...]
                perm_i = sortperm([keep_qubits..., subsystem_qubits...])
                perm_j = sortperm([keep_qubits..., subsystem_qubits...])
                idx_i = sum(full_i[perm_i] .* (2 .^ (0:N-1))) + 1
                idx_j = sum(full_j[perm_j] .* (2 .^ (0:N-1))) + 1
                sum_val += full_density[idx_i, idx_j]
            end
            reduced_matrix[i, j] = sum_val
        end
    end
    
    return reduced_matrix
end

# Function to calculate von Neumann entropy
function von_neumann_entropy(density_matrix::Matrix{Complex{Float64}})
    eigenvalues = eigvals(density_matrix)
    return -sum(λ -> λ > 0 ? λ * log2(λ) : 0, eigenvalues)
end

# Function to calculate entanglement for all bipartitions
function calculate_entanglement(state_vector::Vector{Complex{Float64}})
    N = Int(log2(length(state_vector)))
    entropies = Float64[]
    
    for k in 1:div(N, 2)
        for subsystem in combinations(1:N, k)
            reduced_matrix = reduced_density_matrix(state_vector, collect(subsystem))
            push!(entropies, von_neumann_entropy(reduced_matrix))
        end
    end
    
    return entropies
end

# Function to generate a random quantum state
function random_quantum_state(N::Int)
    state = randn(Complex{Float64}, 2^N)
    return state / norm(state)
end

# Main function to generate states and calculate entanglement
function main(num_states::Int, num_qubits::Int)
    for i in 1:num_states
        state = random_quantum_state(num_qubits)
        entropies = calculate_entanglement(state)
        avg_entropy = mean(entropies)
        max_entropy = maximum(entropies)
        
        println("State $i:")
        println("  Average Entanglement: $avg_entropy")
        println("  Maximum Entanglement: $max_entropy")
        println("  Entropies: $entropies")
        println()
    end
end

# Run the main function
main(5000, 4)  # Generate 5 random states of 4 qubits each



using LinearAlgebra

# Function to generate computational basis states
function generate_computational_basis_states(N::Int)
    states = Vector{Vector{ComplexF64}}()
    for i in 0:(2^N - 1)
        state = zeros(ComplexF64, 2^N)
        state[i + 1] = 1  # Julia uses 1-based indexing
        push!(states, state)
    end
    return states
end

# Function to apply Hadamard gate
function apply_hadamard(state::Vector{ComplexF64}, qubit::Int)
    N = Int(log2(length(state)))
    H = 1/√2 * [1 1; 1 -1]
    return kron(I(2^(N-qubit-1)), H, I(2^qubit)) * state
end

# Function to apply CNOT gate
function apply_cnot(state::Vector{ComplexF64}, control::Int, target::Int)
    N = Int(log2(length(state)))
    CNOT = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]
    if control < target
        return kron(I(2^(N-target-1)), CNOT, I(2^control)) * state
    else
        return kron(I(2^(N-control-1)), CNOT, I(2^target)) * state
    end
end

# Function to generate stabilizer states
function generate_stabilizer_states(N::Int)
    states = generate_computational_basis_states(N)
    
    # Apply Hadamard gates
    for i in 0:(N-1)
        new_states = Vector{Vector{ComplexF64}}()
        for state in states
            push!(new_states, state)
            push!(new_states, apply_hadamard(state, i))
        end
        states = new_states
    end
    
    # Apply CNOT gates
    for control in 0:(N-1)
        for target in 0:(N-1)
            if control != target
                new_states = Vector{Vector{ComplexF64}}()
                for state in states
                    push!(new_states, state)
                    push!(new_states, apply_cnot(state, control, target))
                end
                states = new_states
            end
        end
    end
    
    # Remove duplicates and normalize
    unique_states = unique(states)
    return [state / norm(state) for state in unique_states]
end

# Example usage
N = 1  # Number of qubits
stabilizer_states = generate_stabilizer_states(N)
println("Number of stabilizer states for $N qubits: $(length(stabilizer_states))")

# Print first few states
for (i, state) in enumerate(stabilizer_states[1:min(5, end)])
    println("State $i: ", state)
end

VC = generate_computational_basis_states(1)
VS = generate_stabilizer_states(1)
include("Random_Unitaries.jl")
include("Magic.jl")
using .Random_Unitary_Generation
using .Measure_Magic

U = Random_Unitary_Generation.Generate_Regular_Unitary_Circuit(1);
Strings = Measure_Magic.GenerateAllPauliStrings(1)
PauliOperators = Measure_Magic.PauliOperatorList(Strings, 1)
MC1 = 0.5 * Measure_Magic.MeasureMagic(U * VC[1], PauliOperators, 2)
MC2 = 0.5 * Measure_Magic.MeasureMagic(U * VC[2], PauliOperators, 2)
MC1 + MC2
MS1 = 1/6 * Measure_Magic.MeasureMagic(U * VS[1], PauliOperators, 2)
MS2 = 1/6 * Measure_Magic.MeasureMagic(U * VS[2], PauliOperators, 2)
MS3 = 1/6 * Measure_Magic.MeasureMagic(U * VS[3], PauliOperators, 2)
MS4 = 1/6 * Measure_Magic.MeasureMagic(U * VS[4], PauliOperators, 2)
MS5 = 1/6 * Measure_Magic.MeasureMagic(U * VS[5], PauliOperators, 2)
MS6 = 1/6 * Measure_Magic.MeasureMagic(U * VS[7], PauliOperators, 2)
MS1 + MS2 + MS3 + MS4 + MS5 + MS6
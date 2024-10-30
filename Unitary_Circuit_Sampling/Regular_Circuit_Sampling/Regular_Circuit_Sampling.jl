# Sampling of N qubit random unitary circuits.
# The only parameter here is the size of the circuit, that is 2^N by 2^N
# N - Number of qubits
# 
# The starting state is optional, but due to the CUE matrices used any state suffices.
# |Ψ⟩ = |0⟩
# U|Ψ⟩ = U|0⟩ = |Haar State⟩

using Random

include("Random_Unitaries.jl")
using .Random_Unitary_Generation

# Setting the seed for the random number generation
Seed = 1
Random.seed!(Seed)

# Sampling parameters
No_Samples = 2^20

for i in 1:No_Samples





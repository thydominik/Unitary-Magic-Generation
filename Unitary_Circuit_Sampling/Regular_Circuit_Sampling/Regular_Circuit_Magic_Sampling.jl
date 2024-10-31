# Sampling of N qubit random unitary circuits.
# The only parameter here is the size of the circuit, that is 2^N by 2^N
# N - Number of qubits
# 
# The starting state is optional, but due to the CUE matrices used any state suffices.
# |Ψ⟩ = ⊗_i^N |+>
# U|Ψ⟩ = U|0⟩ = |Haar State⟩

using Random
using ProgressBars

include("Random_Unitaries.jl")
include("Magic.jl")
using .Random_Unitary_Generation
using .Measure_Magic

# Setting the seed for the random number generation
Seed = 1
Random.seed!(Seed)

# Sampling parameters
No_Samples = 2^20

No_Qubits = 2

Psi_0 = 1/sqrt(2^No_Qubits) * ones(2^No_Qubits);


Strings = Measure_Magic.GenerateAllPauliStrings(No_Qubits)
PauliOperators = Measure_Magic.PauliOperatorList(Strings, No_Qubits)

Magic = Vector()

for i in ProgressBar(1:No_Samples)
    t1 = time();
    
    U = Random_Unitary_Generation.Generate_Regular_Unitary_Circuit(No_Qubits);
    State = U * Psi_0
    push!(Magic, Measure_Magic.MeasureMagic(State, PauliOperators, 2))
    # println("Depth = ", D, " sample: ", i, " Time: ", time() - t1)
end



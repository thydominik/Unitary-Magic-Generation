# Sampling of N qubit random unitary Brick Wall circuits.
# There are two paramters to a BW circuit
# D - Depth of the circuit
# N - Number of qubits
# 
# The starting state is irrelevant, but due to the CUE matrices used, any state suffices, I've chosen a 0 magic state.
# |Ψ⟩ = ⊗_i^N |+>
# U|Ψ⟩ = |Haar State⟩

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

# Set the number of qubits first
No_Qubits = 2

# generate the first/Initial state
Psi_0 = 1/sqrt(2^No_Qubits) * ones(2^No_Qubits);

# Generate all possible Pauli strings (as char):
Strings = Measure_Magic.GenerateAllPauliStrings(No_Qubits)
# Predetermine the pauli operators and save them, these will be used at every magic calculation
PauliOperators = Measure_Magic.PauliOperatorList(Strings, No_Qubits)

MaxDepth = No_Qubits

Magic = Vector{Float32}()
for D in ProgressBar(1:MaxDepth)
    for i in ProgressBar(1:No_Samples)
        t1 = time();
        
        U = Random_Unitary_Generation.Generate_BW_Unitary_Circuit(No_Qubits, D);
        State = U * Psi_0
        push!(Magic, Measure_Magic.MeasureMagic(State, PauliOperators, 2))
        # println("Depth = ", D, " sample: ", i, " Time: ", time() - t1)
    end

end


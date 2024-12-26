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
using JLD2
using Base.Threads

current_dir = @__DIR__

# Unitary matrices
filepath = joinpath(current_dir, "..", "Modules", "Magic.jl")
include(filepath)
# Magic
filepath = joinpath(current_dir, "..", "Modules", "Random_Unitaries.jl")
include(filepath)
#include("C:\\Dominik\\PhD\\Projects\\Unitary-Magic-Generation\\Modules\\Magic.jl")

using .Measure_Entanglement
using .Measure_Magic
using .Random_Unitary_Generation

function BrickWall_Circuit_Magic_sampling(No_Qubits::Int, Depth::Int, No_Samples::Int, Seed::Int)
    # Setting the seed for the random number generation
    Random.seed!(Seed)

    # Sampling parameters
    No_Samples = No_Samples

    # Set the number of qubits first
    No_Qubits = No_Qubits

    Depth = Depth
 
    Psi_0 = 1/sqrt(2^No_Qubits) * ones(2^No_Qubits);

    Strings = Measure_Magic.GenerateAllPauliStrings(No_Qubits)
    PauliOperators = Measure_Magic.PauliOperatorList(Strings, No_Qubits)

    Magic           = Vector{Float64}()
    Entanglement    = Vector{Vector{Float64}}()
    SubSystems      = Vector{Vector{Vector{Int}}}()

    for i in ProgressBar(1:No_Samples)    
        U       = Random_Unitary_Generation.Generate_BW_Unitary_Circuit(No_Qubits, Depth);
        State   = U * Psi_0
        
        SVN, SubSys = Measure_Entanglement.calculate_entanglement(State, subsystems = "Middle")
        push!(Magic, Measure_Magic.MeasureMagic_Pure(State, PauliOperators, 2))
    end

    fname = "BWUnitaryCircuitMagicSampled_N_$(No_Qubits)_D_$(Depth)_Samples_$(No_Samples)_Seed_$(Seed).jld2"
    @save fname Magic Psi_0 No_Samples No_Qubits Depth Seed
end


N = 5
for D in [9, 10, 15]
    Partitions = 1
    div = Int(log2(Partitions))
    for s in 1:Partitions
        BrickWall_Circuit_Magic_sampling(N, D, 2^(20-div), s)
    end
end



# Sampling of N qubit random unitary circuits.
# The only parameter here is the size of the circuit, that is 2^N by 2^N
# N - Number of qubits
# 
# The starting state is optional, but due to the CUE matrices used any state suffices.
# |Ψ⟩ = ⊗_i^N |+>
# U|Ψ⟩ = U|0⟩ = |Haar State⟩

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
# Measure_Entanglement
filepath = joinpath(current_dir, "..", "Modules", "Entanglement.jl")
include(filepath)

using .Measure_Entanglement
using .Measure_Magic
using .Random_Unitary_Generation

function Regular_Circuit_Magic_sampling(No_Qubits::Int, No_Samples::Int, Seed::Int)
    # Setting the seed for the random number generation
    Random.seed!(Seed)

    # Sampling parameters
    No_Samples = No_Samples
    # Set the number of qubits first
    No_Qubits = No_Qubits

    Psi_0 = 1/sqrt(2^No_Qubits) * ones(2^No_Qubits);

    Strings = Measure_Magic.GenerateAllPauliStrings(No_Qubits)
    PauliOperators = Measure_Magic.PauliOperatorList(Strings, No_Qubits)

    Magic           = Vector{Float64}()
    Entanglement    = Vector{Vector{Float64}}()
    SubSystems      = Vector{Vector{Vector{Int}}}()

    for i in ProgressBar(1:No_Samples)    
        U = Random_Unitary_Generation.Generate_Regular_Unitary_Circuit(No_Qubits);
        State = U * Psi_0
        SVN, SubSys = Measure_Entanglement.calculate_entanglement(State, subsystems = "Middle")
            # println(typeof(SubSys))
            # println("iteration : ", i)
            # println("Subsystem: ", SubSys)
            # println("Entropies: ", SVN)
        #push!(Entanglement, SVN)
        #push!(SubSystems, SubSys)
    
        push!(Magic, Measure_Magic.MeasureMagic_Pure(State, PauliOperators)[1])
            # println("Depth = ", D, " sample: ", i, " Time: ", time() - t1)
    end

    fname = "RegularUnitaryCircuitMagicSampled_N_$(No_Qubits)_Samples_$(No_Samples)_Seed_$(Seed).jld2"
    #matwrite(fname, Dict("Magic" => Magic)) 
    @save fname Psi_0 No_Samples No_Qubits Seed Entanglement SubSystems
    return Magic, Entanglement
end

N = 7
Partitions = 1
div = Int(log2(Partitions))
M, S = Regular_Circuit_Magic_sampling(N, 2^(20), 2)
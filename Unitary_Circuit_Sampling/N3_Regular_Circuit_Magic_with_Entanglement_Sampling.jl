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
using LinearAlgebra
using Plots
using StatsBase

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

function NQubit_Regular_Circuit_Magic_w_Entanglement_sampling(No_Qubits::Int, No_Samples::Int, Seed::Int)
    # Setting the seed for the random number generation
    Random.seed!(Seed)

    # Sampling parameters
    No_Samples = No_Samples
    # Set the number of qubits first
    No_Qubits = No_Qubits

    Psi_0 = 1/sqrt(2^No_Qubits) * ones(2^No_Qubits);

    PauliOperators = Dict();
    for pauliIndex in 1:No_Qubits
        PauliOperators[pauliIndex] = Measure_Magic.PauliOperatorList(Measure_Magic.GenerateAllPauliStrings(pauliIndex), pauliIndex)
    end

    Magic           = Vector{Float64}()
    SMagic          = Vector{Vector{Float64}}()
    Entanglement    = Vector{Vector{Float64}}()
    SubSystems      = Vector{Vector{Vector{Int}}}()
    Purity          = Vector{Vector{Float64}}()
    for i in ProgressBar(1:No_Samples)    
        U = Random_Unitary_Generation.Generate_Regular_Unitary_Circuit(No_Qubits);
        State = U * Psi_0

        Entropies   = Float64[]
        SubSys      = Vector{Int}[]
        SubMagic    = Vector{Float64}()
        Pure        = Vector{Float64}()
        for q in 1:(No_Qubits - 1)
            subsystem = 1:q
            reduced_matrix_left,    keep_qubits = Measure_Entanglement.reduced_density_matrix(State, collect(subsystem))
            subsystem = (q + 1):No_Qubits
            reduced_matrix_right,   keep_qubits = Measure_Entanglement.reduced_density_matrix(State, collect(subsystem))

            push!(Entropies, Measure_Entanglement.von_neumann_entropy(reduced_matrix_left))
            push!(SubSys, keep_qubits)
            push!(SubMagic, Measure_Magic.MeasureMagic_Mixed(reduced_matrix_left,   PauliOperators[No_Qubits - q], 2))
            push!(SubMagic, Measure_Magic.MeasureMagic_Mixed(reduced_matrix_right,  PauliOperators[q], 2))
            push!(Pure, real(sum(diag(reduced_matrix_left^2))))
            push!(Pure, real(sum(diag(reduced_matrix_right^2))))
        end

        push!(Entanglement, Entropies)
        push!(SubSystems, SubSys)
        push!(Magic, Measure_Magic.MeasureMagic_Pure(State, PauliOperators[No_Qubits], 2))
        push!(SMagic, SubMagic)
        push!(Purity, Pure)
            # println("Depth = ", D, " sample: ", i, " Time: ", time() - t1)
    end

    #fname = "RegularUnitaryCircuitMagicSampled_N_$(No_Qubits)_Samples_$(No_Samples)_Seed_$(Seed).jld2"
    #matwrite(fname, Dict("Magic" => Magic)) 
    #Psi_0 #No_Samples #No_Qubits #Seed #Entanglement #SubSystems
    return Magic, Entanglement, SMagic, SubSystems, Purity
end


M, S, SM, SS, Pu = NQubit_Regular_Circuit_Magic_w_Entanglement_sampling(1, 2^20, 2);
SML = map(first, SM)
SMR = map(last, SM)
S   = map(first, S)
Pu  = map(first, Pu)

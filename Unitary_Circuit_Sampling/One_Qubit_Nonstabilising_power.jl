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
using MAT

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

Seed = 1
# Setting the seed for the random number generation
Random.seed!(Seed)
# Sampling parameters
No_Samples = 2^20
No_Stab = 6
# Set the number of qubits first
No_Qubits = 1

Psi = zeros(Complex, 6, 2);
Psi[1, :] = [0 1];
Psi[2, :] = [1 0];
Psi[3, :] = [1/sqrt(2) 1/sqrt(2)];
Psi[4, :] = [1/sqrt(2) -1/sqrt(2)];
Psi[5, :] = [1/sqrt(2) im/sqrt(2)];
Psi[6, :] = [1/sqrt(2) -im/sqrt(2)];
Strings = Measure_Magic.GenerateAllPauliStrings(No_Qubits)
PauliOperators = Measure_Magic.PauliOperatorList(Strings, No_Qubits)
Magic           = Vector{Float64}()
#Entanglement    = Vector{Vector{Float64}}()
#SubSystems      = Vector{Vector{Vector{Int}}}()

for i in ProgressBar(1:No_Samples) 
    M2 = Vector{Float64}()
    for stabstate in 1:No_Stab
        Psi_0 = Psi[stabstate, :] 
        U = Random_Unitary_Generation.Generate_Regular_Unitary_Circuit(No_Qubits);
        State = U * Psi_0

        push!(M2, Measure_Magic.MeasureMagic(State, PauliOperators, 2))
    end
    push!(Magic, sum(M2)/No_Stab)
end
fname = "NonStabilisingPower_of_Unitary_N_1_Samples_$(No_Samples)_Seed_$(Seed).mat"
matwrite(fname, Dict("Magic" => Magic))



using Plots
using Statistics
plot(sort(Magic))
using StatsBase
histogram(Magic, nbins = 500)
mean(Magic)

skewness(Magic)
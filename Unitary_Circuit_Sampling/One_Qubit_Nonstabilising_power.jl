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
No_Samples = 2^22
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
Magic           = zeros(No_Samples, No_Stab)


for i in ProgressBar(1:No_Samples) 
    #M2 = Vector{Float64}()
    U = Random_Unitary_Generation.Generate_Regular_Unitary_Circuit(No_Qubits);
    for stabstate in 1:No_Stab
        Psi_0 = Psi[stabstate, :] 
        
        State = U * Psi_0

        #push!(M2, Measure_Magic.MeasureMagic_Pure(State, PauliOperators, 2)[1])
        Magic[i, stabstate] = Measure_Magic.MeasureMagic_Pure(State, PauliOperators, 2)[1]
    end
    #push!(Magic, mean(M2))
    #empty!(M2)
end
#fname = "NonStabilisingPower_of_Unitary_N_1_Samples_$(No_Samples)_Seed_$(Seed).mat"
#matwrite(fname, Dict("Magic" => Magic))

fname = "UnitaryMagicSampling_AllStabs_N_1_Samples_$(No_Samples).jld2"
@save fname Magic

using Plots
plot()
histogram!(Magic[:, 1], bins=range(0, log2(3/2), 1000), normalize=:pdf)

NSP = zeros(No_Samples * 2)
for i in 1:No_Samples*2
    NSP[i] = sum(Magic[i, :]) / No_Stab
end

histogram!(NSP, bins=range(0, log2(3/2), 1000), normalize=:pdf)

plot!(xlims=(0, log2(3/2)), ylims=(0, 7))
plot(size=(400, 400), Magic[1:100, 3], Magic[1:100, 6], xlims=(0, log2(3/2)), ylims=(0, log2(3/2)))


vline!([mean(Magic), mean(NSP)])
mean(Magic)
mean(NSP)
using Statistics
plot(sort(Magic)[1:10:end])
using StatsBase
histogram(Magic, nbins = 500, normalize=:pdf)
mean(Magic)

skewness(Magic)




using LaTeXStrings
Magic = load("UnitaryMagicSampling_AllStabs_N_1_Samples_4194304.jld2")["Magic"]
NSP = zeros(No_Samples * 2)
for i in 1:No_Samples*2
    NSP[i] = sum(Magic[i, :]) / 6
end
plot()
histogram!(Magic[:], bins=range(0, log2(3/2), 1000), normalize=:pdf, label=L"\rho(M_2)")
histogram!(NSP, bins=range(0, log2(3/2), 1000), normalize=:pdf, label=L"\rho(\mathcal{M}_2)")
plot!(xlims=(0, log2(3/2)),
    ylims=(0, 7),
    xlabel=L"M_2",
    ylabel=L"\rho(M_2), \rho(\mathcal{M}_2)")
vline!([mean(Magic)], color=:black, label = "Average")

savefig("UnitaryMagicSampling_Histograms.pdf")


using Statistics
using Plots
using StatsBase
using JLD2

No_Qubits = 3
Depth = 2
fname = "BWUnitaryCircuitMagicSampled_N_$(No_Qubits)_D_$(Depth)_Samples_1048576_Seed_1.jld2"
fname = "C:/Dominik/PhD/Projects/Unitary-Magic-Generation/Data/RegularUnitaryCircuits/RegularUnitaryCircuitMagicSampled_N_$(No_Qubits)_Samples_1048576_Seed_1.jld2"
#fname = "C:/Dominik/PhD/Projects/Unitary-Magic-Generation/Data/RegularUnitaryCircuits/RegularUnitaryCircuitMagicSampled_N_$(No_Qubits)_Samples_1024_Seed_1_no_Ent.jld2"

data    = JLD2.load(fname)
M2      = data["Magic"]
plot(sort(M2)/(log((2^No_Qubits + 1)/2)))

a = quantile(M2, [0.25, 0.5, 0.75])
edges = 1000
h = fit(Histogram, M2, edges=[0:0.01:10])

bins = Vector(h.edges[1])
counts = h.weights

plot(bins[1:length(bins)-1], counts)
histogram!(M2)
vline!(a)
vline!([mean(M2)])
vline!([median(M2)])

mean(M2)/log((2^No_Qubits + 1)/2)

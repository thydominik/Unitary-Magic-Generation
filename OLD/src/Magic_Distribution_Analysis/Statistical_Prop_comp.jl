using Statistics
using Plots
using PlotThemes
using LaTeXStrings
using StatsBase
using JLD2


Statdata = Dict()
for keys in range(1,5)
    Statdata[keys] = JLD2.load("N$(keys)_Statistics.jld2")
end

data = Dict()
for i in range(1, 5)
    println(i)
    dataPath = "D:\\Data\\Random_Unitary_Magic_Generation\\RegularUnitaryCircuitMagicSampled_N_$(i)_Samples_1048576_Seed_1.jld2"
    data[i] = JLD2.load(dataPath)
end

#  Distribution comparison
No_Data = 5
M2 = Dict()
for i in range(1, 5)
    M2[i] = round.(data[i]["Magic"] ./ log((2^i + 1)/2),    digits=14)
    #M2[i] = round.(data[i]["Magic"],    digits=14)
end

h = Dict()
for i in range(1, No_Data)
    bin_edges = range(0, 1, length = 1001)
    h[i] = fit(Histogram, M2[i], bin_edges)
    h[i] = StatsBase.normalize(h[i], mode=:pdf)
end

counts = Dict()
for i in range(1, No_Data)
    counts[i] = h[i].weights
end

p = plot(dpi=400)

for keys in range(1, 5)
    #plot!(range(0, log((2^keys + 1)/2), length = 1000), counts[keys], lw = 2, label="N = $(keys)", legend=:topright)
    plot!(bin_edges[1:end-1], counts[keys], lw = 2, label="N = $(keys)", legend=:topleft)
end
display(p)
theme(:mute::Symbol;)

title!(L"Magic Distribution $N = 1-5$",   titlefontsize=20)
xlabel!(L"$\tilde{M}_2$",               labelfontsize=20)
ylabel!(L"$\varrho(\tilde{M}_2)$",      labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)
vline!([mean(M2[1]) mean(M2[2]) mean(M2[3]) mean(M2[4]) mean(M2[5])], color="black", label = "averages")
savefig(p, "n_1_5_normalised.pdf")
savefig(p, "n_1_5_normalised.png")

# Mean Magic ===============================================================
MeanMagic = Vector()
for keys in range(1, 5)
    if keys == 1 || keys == 2
        push!(MeanMagic, Statdata[keys]["Mean"][end])
    else
        push!(MeanMagic, maximum(Statdata[keys]["Mean"]))
    end
end
N = [1, 2, 3, 4 ,5]
p = plot(dpi=400)
scatter!(N, [mean(M2[1]), mean(M2[2]), mean(M2[3]), mean(M2[4]), mean(M2[5])], label="Mean")

# Median Magic ===============================================================

scatter!(N, [median(M2[1]), median(M2[2]), median(M2[3]), median(M2[4]), median(M2[5])], label="Median")
plot!(framestyle=:box)
plot!(legendfontsize=10)
#title!(L"Magic Distribution $N = 1-5$",   titlefontsize=20)
xlabel!(L"N",               labelfontsize=20)
ylabel!("Statstical Vars",      labelfontsize=20)
savefig(p, "n_1_5_mean_median.pdf")
savefig(p, "n_1_5_mean_median.png")

p = plot(dpi=400)
scatter!(N, [skewness(M2[1]), skewness(M2[2]), skewness(M2[3]), skewness(M2[4]), skewness(M2[5])], label="Skewness")
scatter!(N, [kurtosis(M2[1]), kurtosis(M2[2]), kurtosis(M2[3]), kurtosis(M2[4]), kurtosis(M2[5])], label="Kurtosis")
plot!(framestyle=:box)
plot!(legendfontsize=10)
#title!(L"Magic Distribution $N = 1-5$",   titlefontsize=20)
xlabel!(L"N",               labelfontsize=20)
ylabel!("Statstical Vars",      labelfontsize=20)
savefig(p, "n_1_5_skew_kurt.pdf")
savefig(p, "n_1_5_skew_kurt.png")
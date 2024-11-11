using Statistics
using Plots
using PlotThemes
using LaTeXStrings
using StatsBase
using JLD2

No_Qubits = 5
No_Depth = No_Qubits + 1
No_Data = No_Depth + 1
# Example:
#plot(x, y, title="Customized Plot", xlabel="X-axis", ylabel="Y-axis", lw=3, color=:red, marker=:circle, legend=:topright)

# Loading data -------------------------------------------------------------------------------------------------------------------------------------------
data = Dict()
dataPath = "D:\\Data\\Random_Unitary_Magic_Generation\\RegularUnitaryCircuitMagicSampled_N_$(No_Qubits)_Samples_1048576_Seed_1.jld2"
data[1]    = JLD2.load(dataPath)
for i in range(1, No_Depth)
    dataPath = "D:\\Data\\Random_Unitary_Magic_Generation\\BWUnitaryCircuitMagicSampled_N_$(No_Qubits)_D_$(i)_Samples_1048576_Seed_1.jld2"
    data[i + 1] = JLD2.load(dataPath)
end

# Manipulation of data -------------------------------------------------------------------------------------------------------------------------------------------

M2 = Dict()
for i in range(1, No_Data)
    M2[i] = round.(data[i]["Magic"] ./ log((2^No_Qubits + 1)/2),    digits=14)
end


# Quuantiles
Quantiles = Dict()
for i in range(1, No_Data)
    Quantiles[i] = quantile(M2[i], [0.25, 0.5, 0.75])
end

# Create histogram for the normalised magic data
bin_edges = range(0, 1, length = 1001)

h = Dict()
for i in range(1, No_Data)
    h[i] = fit(Histogram, M2[i], bin_edges)
    h[i] = StatsBase.normalize(h[i], mode=:pdf)
end

counts = Dict()
for i in range(1, No_Data)
    counts[i] = h[i].weights
end

# Check normalisation
Integral = Dict()
for i in range(1, No_Data)
    Integral[i] = 0
end

for i in range(1, 5), j in 1:(length(bin_edges) - 1)
    Integral[i] += (bins[j + 1] - bins[j]) * counts[i][j]
end
println(Integral)


# plotting -------------------------------------------------------------------------------------------------------------------------------------------
p = plot(dpi=400)
for keys in range(1, No_Data)
    if keys == 1
        plot!(bin_edges[1:length(bin_edges)-1], (counts[keys]), lw=2, label="Regular circuit", legend=:topleft)
    else
        plot!(bin_edges[1:length(bin_edges)-1], (counts[keys]), lw=2, label="D = $(keys - 1)", legend=:topleft)
    end
end
display(p)
theme(:mute::Symbol;)

title!(L"Magic Distribution $N = 5$ and $D = 1-6$",   titlefontsize=20)
xlabel!(L"$\tilde{M}_2$",               labelfontsize=20)
ylabel!(L"$\varrho(\tilde{M}_2)$",      labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)

savefig(p, "N$(No_Qubits)_Depth_vs_reg.pdf")
savefig(p, "N$(No_Qubits)_Depth_vs_reg.png")

# plotting -------------------------------------------------------------------------------------------------------------------------------------------
p = plot(dpi=400)
plot!(bin_edges[1:length(bin_edges)-1], (counts[1]), lw=2, label="Magic Distribution", legend=:topleft)
vline!([Quantiles[1][1]], color=:red,      label="25% = $(round(Quantiles[1][1], digits=5))")
vline!([Quantiles[1][2]], color=:green,    label="50% = $(round(Quantiles[1][2], digits=5))")
vline!([Quantiles[1][3]], color=:blue,     label="75% = $(round(Quantiles[1][3], digits=5))")
vline!([mean(M2[1])],      linestyle=:dash, label="Average Magic = $(round(mean(M2[1]),   digits=4))")
vline!([median(M2[1])],    linestyle=:dash, label="Median Magic = $(round(median(M2[1]),  digits=4))")

theme(:mute::Symbol;)

title!(L"Magic Distribution $N = 4$",   titlefontsize=20)
xlabel!(L"$\tilde{M}_2$",               labelfontsize=20)
ylabel!(L"$\varrho(\tilde{M}_2)$",      labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)

annotate!((bin_edges[argmax(counts[1])], 2.5, text("max(P(M)) = $(bin_edges[argmax(counts[1])])"), :center))

savefig(p, "N$(No_Qubits)_RegularCircuit_MagicDistribution.pdf")
savefig(p, "N$(No_Qubits)_RegularCircuit_MagicDistribution.png")

# plotting -------------------------------------------------------------------------------------------------------------------------------------------

p = plot(dpi=400)
plot!(bin_edges[1:length(bin_edges)-1], (counts[2]), lw=6, label="Magic Distribution N = 5, D = 1", legend=:bottomright)
h2 = fit(Histogram, M2_2["Magic"] ./ log((2^No_Qubits + 1)/2), bin_edges)
h2 = StatsBase.normalize(h2, mode=:pdf)

plot!(bin_edges[1:length(bin_edges)-1], h2.weights, lw=1, label="Magic Distribution N = 2", legend=:bottomright)

title!(L"Magic Distribution $N = 3, D = 1$ & $N = 2 $",   titlefontsize=20)
xlabel!(L"$\tilde{M}_2$",               labelfontsize=20)
ylabel!(L"$\varrho(\tilde{M}_2)$",      labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)

savefig(p, "N3_D1_and_N2_compariosn.pdf")
savefig(p, "N3_D1_and_N2_compariosn.png")


# Statistical tests:
samples     = zeros(No_Data, 20)
Mean        = zeros(No_Data, 20)
Median      = zeros(No_Data, 20)
Skewness    = zeros(No_Data, 20)
Kurtosis    = zeros(No_Data, 20)
Variance    = zeros(No_Data, 20)
k = 1
for keys in range(1, No_Data), i in 1:1:20
    println(keys)
    samples[keys, i]    = 2^i
    Mean[keys, i]       = mean(M2[keys][1:Int(samples[keys, i])])
    Median[keys, i]     = median(M2[keys][1:Int(samples[keys, i])])
    Skewness[keys, i]   = skewness(M2[keys][1:Int(samples[keys, i])])   
    Kurtosis[keys, i]   = kurtosis(M2[keys][1:Int(samples[keys, i])])
    Variance[keys, i]   = var(M2[keys][1:Int(samples[keys, i])])
end


@save "N$(No_Qubits)_Statistis.jld2" counts bin_edges samples Mean Median Skewness kurtosis Variance

# plotting -------------------------------------------------------------------------------------------------------------------------------------------

p = plot(dpi=400)
for keys in range(2, No_Data)
    scatter!(log2.(samples[keys, :]), Mean[keys, :], label="Mean, D = $(keys - 1)")
end
scatter!(log2.(samples[1, :]), Mean[1, :], label="Mean - reg. Circ")
display(p)
theme(:mute::Symbol;)

title!(L"Mean $N = 5$", titlefontsize=20)
xlabel!("Sample size",      labelfontsize=20)
ylabel!(L"$\tilde{M}_2$",           labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)

savefig(p, "N$(No_Qubits)_sample_vs_mean.pdf")
savefig(p, "N$(No_Qubits)_sample_vs_mean.png")

# plotting -------------------------------------------------------------------------------------------------------------------------------------------

p = plot()
for keys in range(2, No_Data)
    scatter!(log2.(samples[keys, :]), Median[keys, :], label="Median, D = $(keys - 1)")
end
scatter!(log2.(samples[1, :]), Mean[1, :], label="Median - reg. Circ")
display(p)
theme(:mute::Symbol;)

title!(L"Median $N = 5$", titlefontsize=20)
xlabel!("Sample size",      labelfontsize=20)
ylabel!(L"$\tilde{M}_2$",           labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)

savefig(p, "N$(No_Qubits)_sample_vs_median.pdf")
savefig(p, "N$(No_Qubits)_sample_vs_median.png")

# plotting -------------------------------------------------------------------------------------------------------------------------------------------

p = plot(dpi=400)
for keys in range(2, No_Data)
    scatter!(log2.(samples[keys, :]), Variance[keys, :], label="Var, D = $(keys - 1)")
end
scatter!(log2.(samples[1, :]), Variance[1, :], label="Var - reg. Circ")
display(p)
theme(:mute::Symbol;)

title!(L"Variance $N = 5$", titlefontsize=20)
xlabel!("Sample size",      labelfontsize=20)
ylabel!(L"$\tilde{M}_2$",           labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)

savefig(p, "N$(No_Qubits)_sample_vs_var.pdf")
savefig(p, "N$(No_Qubits)_sample_vs_var.png")

# plotting -------------------------------------------------------------------------------------------------------------------------------------------

p = plot(dpi=400)
for keys in range(2, No_Data)
    scatter!(log2.(samples[keys, :]), Skewness[keys, :], label="Skewness, D = $(keys - 1)")
end
scatter!(log2.(samples[1, :]), Skewness[1, :], label="Skewness- reg. Circ")
display(p)
theme(:mute::Symbol;)

title!(L"Skewness $N = 5$", titlefontsize=20)
xlabel!("Sample size",      labelfontsize=20)
ylabel!(L"$\tilde{M}_2$",           labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)

savefig(p, "N$(No_Qubits)_sample_vs_skew.pdf")
savefig(p, "N$(No_Qubits)_sample_vs_skew.png")

# plotting -------------------------------------------------------------------------------------------------------------------------------------------

p = plot(dpi=400)
for keys in range(2, No_Data)
    scatter!(log2.(samples[keys, :]), Kurtosis[keys, :], label="Kurtosis, D = $(keys - 1)")
end
scatter!(log2.(samples[1, :]), Kurtosis[1, :], label="Kurtosis - reg. Circ", legend=:bottomright)
display(p)
theme(:mute::Symbol;)

title!(L"Kurtosis $N = 5$", titlefontsize=20)
xlabel!("Sample size",      labelfontsize=20)
ylabel!(L"$\tilde{M}_2$",           labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)

savefig(p, "N$(No_Qubits)_sample_vs_Kurtosis.pdf")
savefig(p, "N$(No_Qubits)_sample_vs_Kurtosis.png")


# Entanglement VS Magic ============================================================
using KernelDensity
using StatsPlots

SVN = Vector{Float64}()
for i in 1:(2^20)
    push!(SVN, round(mean(data[1]["Entanglement"][i][:]), digits=14))
end

p = plot(dpi=400, legend=:topright)
# scatter!(M2[1:10000], SVN[1:10000])

k = kde((M2[1], SVN))
p = contourf(k, c = :vik, linewidth = 1)

theme(:mute::Symbol;)

title!(L"Kernel Density: $mean(S_{VN})$ vs $M_2$ $N = 4$", titlefontsize=20)
xlabel!(L"$\tilde{M}_2$",      labelfontsize=20)
ylabel!(L"$S_{VN}$",           labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)

plot!(xscale=:lin)
#plot!(xlims=[0, 1])
#plot!(ylims=[0, log2(2.5)])
savefig(p, "N$(No_Qubits)_SVN_vs_M2_mean.pdf")
savefig(p, "N$(No_Qubits)_SVN_vs_M2_mean.png")
using Statistics
using Plots
using PlotThemes
using LaTeXStrings
using StatsBase
using JLD2


# Example:
#plot(x, y, title="Customized Plot", xlabel="X-axis", ylabel="Y-axis", lw=3, color=:red, marker=:circle, legend=:topright)

# Loading data -------------------------------------------------------------------------------------------------------------------------------------------
#dataPath = "D:\\Data\\Random_Unitary_Magic_Generation\\RegularUnitaryCircuitMagicSampled_N_1_Samples_1048576_Seed_1.jld2"
dataPath    = "RegularUnitaryCircuitMagicSampled_N_1_Samples_4194304_Seed_1.jld2"
data        = JLD2.load(dataPath)

# Manipulation of data -------------------------------------------------------------------------------------------------------------------------------------------
    # Normalise data:
M2 = round.(data["Magic"] ./ log(3/2), digits=14)
    # calculate the quantiles of the data
a = quantile(M2, [0.25, 0.5, 0.75])
    # Create histogram for the normalised magic data
h = fit(Histogram, M2, nbins=1000)
    # Normalise magic data
h       = StatsBase.normalize(h, mode=:pdf)
bins    = Vector(h.edges[1])
counts  = h.weights

# Check normalisation
Integral = 0
for i in 1:length(counts)
    Integral += (bins[i + 1] - bins[i]) * counts[i]
end
println(round(Integral, digits=10))

# plotting -------------------------------------------------------------------------------------------------------------------------------------------
p = plot(bins[1:length(bins)-1], (counts), lw=2, label="Magic distribution", legend=:topleft, dpi=400)
vline!([a[1]], color=:red,      label="25% = $(round(a[1], digits=5))")
vline!([a[2]], color=:green,    label="50% = $(round(a[2], digits=5))")
vline!([a[3]], color=:blue,     label="75% = $(round(a[3], digits=5))")
vline!([mean(M2)],      linestyle=:dash, label="Average Magic = $(round(mean(M2),   digits=4))")
vline!([median(M2)],    linestyle=:dash, label="Median Magic = $(round(median(M2),  digits=4))")

theme(:mute::Symbol;)

title!(L"Magic Distribution $N = 1$",   titlefontsize=20)
xlabel!(L"$\tilde{M}_2$",               labelfontsize=20)
ylabel!(L"$\varrho(\tilde{M}_2)$",      labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)

annotate!((bins[argmax(counts)], 2.5, text("max(P(M)) = $(bins[argmax(counts)])"), :center))

savefig(p, "N1_RegularCircuit_MagicDistribution.pdf")
savefig(p, "N1_RegularCircuit_MagicDistribution.png")





# Statistical tests:
samples     = Vector()
Mean        = Vector()
Median      = Vector()
Skewness    = Vector()
Kurtosis    = Vector()
Variance    = Vector()
k = 1
for i in range(1,22)
    M2      = round.(data["Magic"][1:2^i] ./ log(3/2), digits=14)
        # calculate the quantiles of the data
    a = quantile(M2, [0.25, 0.5, 0.75])
        # Create histogram for the normalised magic data
    h = fit(Histogram, M2, nbins=1000)
        # Normalise magic data
    h       = StatsBase.normalize(h, mode=:pdf)
    bins    = Vector(h.edges[1])
    counts  = h.weights

    println(2^i)
    push!(samples, 2^i)
    push!(Mean, mean(M2))
    push!(Median, median(M2))
    push!(Skewness, skewness(M2))
    push!(Kurtosis, kurtosis(M2))
    push!(Variance, var(M2))
    k += 1
end

@save "N1_Statistics.jld2" counts bins samples Mean Median Skewness kurtosis Variance

# Mean ------------------------------------------------------------------------------------------------------------------------------------------

p = scatter(log2.(samples),    Mean,       label="Average")
theme(:mute::Symbol;)

title!(L"Mean $N = 1$", titlefontsize=20)
xlabel!(L"$\log_2$ - Sample size",      labelfontsize=20)
ylabel!(L"$\langle \tilde{M}_2 \rangle$",           labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)
plot!(ylims=[0, 1])
hline!([Mean[end]], label="")
#plot!(xscale=:log)
savefig(p, "N1_sample_vs_mean.pdf")
savefig(p, "N1_sample_vs_mean.png")

# Median------------------------------------------------------------------------------------------------------------------------------------------

p = scatter(log2.(samples),    Median,     label="Median")
theme(:mute::Symbol;)

title!(L"Median $N = 1$", titlefontsize=20)
xlabel!(L"$\log_2$ - Sample size",      labelfontsize=20)
ylabel!(L"$Median(\tilde{M}_2)$",           labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)
plot!(ylims=[0, 1])
hline!([Median[end]], label="")

savefig(p, "N1_sample_vs_median.pdf")
savefig(p, "N1_sample_vs_median.png")

# Skewness------------------------------------------------------------------------------------------------------------------------------------------

p = scatter(log2.(samples),    Skewness,   label="Skewness")
theme(:mute::Symbol;)

title!(L"Skewness $N = 1$", titlefontsize=20)
xlabel!(L"$\log_2$ - Sample size",      labelfontsize=20)
ylabel!(L"$Sk(\tilde{M}_2)$",           labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)
plot!(ylims=[-0.5, 0.5])

savefig(p, "N1_sample_vs_skewness.pdf")
savefig(p, "N1_sample_vs_skewness.png")

# Kurtosis------------------------------------------------------------------------------------------------------------------------------------------

p = scatter(log2.(samples),    Kurtosis,   label="Kurtosis")
theme(:mute::Symbol;)

title!(L"Kurtosis $N = 1$", titlefontsize=20)
xlabel!(L"$\log_2$ - Sample size",      labelfontsize=20)
ylabel!(L"$Krt(\tilde{M}_2)$",           labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)
plot!(ylims=[-2.5, 0])

savefig(p, "N1_sample_vs_kurtosis.pdf")
savefig(p, "N1_sample_vs_kurtosis.png")

# Variance ------------------------------------------------------------------------------------------------------------------------------------------

p = scatter(log2.(samples),    Variance,   label="Kurtosis")
theme(:mute::Symbol;)

title!(L"Variance $N = 1$", titlefontsize=20)
xlabel!(L"$\log_2$ - Sample size",      labelfontsize=20)
ylabel!(L"$Var(\tilde{M}_2)$",           labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)
plot!(ylims=[0, 0.15])

savefig(p, "N1_sample_vs_variance.pdf")
savefig(p, "N1_sample_vs_variance.png")
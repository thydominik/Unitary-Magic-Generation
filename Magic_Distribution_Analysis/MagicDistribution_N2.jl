using Statistics
using Plots
using PlotThemes
using LaTeXStrings
using StatsBase
using JLD2


# Example:
#plot(x, y, title="Customized Plot", xlabel="X-axis", ylabel="Y-axis", lw=3, color=:red, marker=:circle, legend=:topright)

# Loading data -------------------------------------------------------------------------------------------------------------------------------------------
dataPath = "D:\\Data\\Random_Unitary_Magic_Generation\\RegularUnitaryCircuitMagicSampled_N_2_Samples_1048576_Seed_1.jld2"
data    = JLD2.load(dataPath)

# Manipulation of data -------------------------------------------------------------------------------------------------------------------------------------------
    # Normalise data:
M2      = round.(data["Magic"] ./ log(5/2), digits=14)
    # calculate the quantiles of the data
a = quantile(M2, [0.25, 0.5, 0.75])
    # Create histogram for the normalised magic data
    bin_edges = range(0, 1, length = 1001)
h = fit(Histogram, M2, bin_edges)
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

title!(L"Magic Distribution $N = 2$",   titlefontsize=20)
xlabel!(L"$\tilde{M}_2$",               labelfontsize=20)
ylabel!(L"$\varrho(\tilde{M}_2)$",      labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)

annotate!((bins[argmax(counts)], 2.5, text("max(P(M)) = $(bins[argmax(counts)])"), :center))

savefig(p, "N2_RegularCircuit_MagicDistribution.pdf")
savefig(p, "N2_RegularCircuit_MagicDistribution.png")

# Statistical tests:
samples     = Vector()
Mean        = Vector()
Variance    = Vector()
Median      = Vector()
Skewness    = Vector()
Kurtosis    = Vector()
k = 1
for i in range(1, 20)
    M2      = round.(data["Magic"][1:2^i] ./ log(5/2), digits=14)
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

@save "N2_Statistics.jld2" counts bins samples Mean Median Skewness kurtosis Variance

# Mean -------------------------------------------------------------------------------

p = scatter(log2.(samples),    Mean,       label="Average", dpi = 400)
theme(:mute::Symbol;)

title!(L"Mean $N = 2$", titlefontsize=20)
xlabel!(L"$\log_2$ - Sample size",      labelfontsize=20)
ylabel!(L"$\tilde{M}_2$",           labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)

savefig(p, "N2_sample_vs_mean.pdf")
savefig(p, "N2_sample_vs_mean.png")

# Mean -------------------------------------------------------------------------------

p = scatter(log2.(samples),    Median,     label="Median", dpi = 400)
theme(:mute::Symbol;)

title!(L"Median $N = 2$", titlefontsize=20)
xlabel!(L"$\log_2$ - Sample size",      labelfontsize=20)
ylabel!(L"$\tilde{M}_2$",           labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)

savefig(p, "N2_sample_vs_median.pdf")
savefig(p, "N2_sample_vs_median.png")

# Mean -------------------------------------------------------------------------------

p = scatter(log2.(samples),    Skewness,   label="Skewness", dpi=400)
theme(:mute::Symbol;)

title!(L"Skewness $N = 2$", titlefontsize=20)
xlabel!(L"$\log_2$ - Sample size",      labelfontsize=20)
ylabel!(L"$\tilde{M}_2$",           labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)

savefig(p, "N2_sample_vs_skewness.pdf")
savefig(p, "N2_sample_vs_skewness.png")

# Mean -------------------------------------------------------------------------------

p = scatter(log2.(samples),    Kurtosis,   label="Kurtosis")
theme(:mute::Symbol;)

title!(L"Kurtosis $N = 2$", titlefontsize=20)
xlabel!(L"$\log_2$ - Sample size",      labelfontsize=20)
ylabel!(L"$\tilde{M}_2$",           labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)

savefig(p, "N2_sample_vs_kurtosis.pdf")
savefig(p, "N2_sample_vs_kurtosis.png")

# Variance -------------------------------------------------------------------------------

p = scatter(log2.(samples),    Variance,   label="Kurtosis")
theme(:mute::Symbol;)

title!(L"Variance $N = 2$", titlefontsize=20)
xlabel!(L"$\log_2$ - Sample size",      labelfontsize=20)
ylabel!(L"$\tilde{M}_2$",           labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)

savefig(p, "N2_sample_vs_variance.pdf")
savefig(p, "N2_sample_vs_variance.png")


# Entanglement VS Magic ============================================================
using KernelDensity
using StatsPlots

M2      = round.(data["Magic"] ./ log(5/2), digits=14)
SVN = Vector{Float64}()
for i in 1:(2^20)
    push!(SVN, round(data["Entanglement"][i][1], digits=14))
end

plot(legend=:topright)
# scatter!(M2[1:10000], SVN[1:10000])

k = kde((M2, SVN))
p = contourf(k, c = :vik, linewidth = 1, dpi=400)

theme(:mute::Symbol;)

title!(L"Kernel Density: $S_{VN}$ vs $M_2$ $N = 2$", titlefontsize=20)
xlabel!(L"$\tilde{M}_2$",      labelfontsize=20)
ylabel!(L"$S_{VN}$",           labelfontsize=20)

plot!(framestyle=:box)
plot!(legendfontsize=10)

plot!(xscale=:lin)
plot!(xlims=[0, 1])
plot!(ylims=[0, log2(2)])
savefig(p, "KDE_N2_SVN_vs_M2.pdf")
savefig(p, "KDE_N2_SVN_vs_M2.png")
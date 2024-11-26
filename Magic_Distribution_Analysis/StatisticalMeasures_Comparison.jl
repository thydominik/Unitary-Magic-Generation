using Statistics
using Plots
using PlotThemes
using LaTeXStrings
using StatsBase
using JLD2
using LinearAlgebra
using ProgressBars
# N = 1-6 data load

data = Dict()
for i in ProgressBar(range(1, 5))
    println(i)
    dataPath = "D:\\Data\\Random_Unitary_Magic_Generation\\N$(i)\\Regular\\RegularUnitaryCircuitMagicSampled_N_$(i)_Samples_1048576_Seed_1.jld2"
    data[i] = JLD2.load(dataPath)
end

dataPath = "D:\\Data\\Random_Unitary_Magic_Generation\\N6\\Regular\\RegularUnitaryCircuitMagicSampled_N_6_Samples_1048576_MultiSeed_w_Ent.jld2"
data[6] = JLD2.load(dataPath)

# Going from log to log2
Magic = Dict()
for i in 1:6
    Magic[i] = -log2.(exp.(-(data[i]["Magic"] .+ log(2^i)))) .-log2(2^i);
end

# Distributions of magic
M = 1000;
h = Dict()
for NIndex in 1:6
    bin_edges = range(0, log((2^NIndex + 1)/2), M)

    if NIndex < 6
        h[NIndex] = fit(Histogram, data[NIndex]["Magic"], bin_edges)
    else
        h[NIndex] = fit(Histogram, data[NIndex]["Magic"], bin_edges)
    end
end

counts = Dict()
for i in range(1, 6)
    counts[i] = h[i].weights
end

# Not scaled, normalised
MostProbable = Vector{Float64}()
p = plot(dpi=400)
for i in 1:6
    bin_edges = range(0, log2((2^i + 1)/2), M-1)
    integral = 0
    for j in 1:length(bin_edges)
        integral += (bin_edges[2] - bin_edges[1]) * counts[i][j]
    end
    push!(MostProbable, bin_edges[argmax(counts[i])])
    plot!(bin_edges, counts[i] ./ integral, label = "N = $(i)")
end
title!("Magic distribution for N = 1-6")
xlabel!(L"$M_2$",   labelfontsize = 20)
ylabel!(L"$\varrho ( M_2 )$",   labelfontsize = 20)

plot!(framestyle=:box)
plot!(legendfontsize=10)

display(p)
savefig(p, "N1_N6_magic_distribution.pdf")
savefig(p, "N1_N6_magic_distribution.png")

# scaled and normalised

p = plot(dpi=400)
for i in 1:6
    bin_edges = range(0, 1, M-1)
    integral = 0
    for j in 1:length(bin_edges)
        integral += (bin_edges[2] - bin_edges[1]) * counts[i][j]
    end

    plot!(bin_edges, counts[i] ./ integral, lw = 2, label = "N = $(i)")
end
title!("Normalised Magic distribution for N = 1-6")
xlabel!(L"$\tilde{M}_2$",   labelfontsize = 20)
ylabel!(L"$\varrho ( \tilde{M}_2 )$",   labelfontsize = 20)

plot!(framestyle=:box)
plot!(legendfontsize=10)

display(p)
savefig(p, "N1_N6_magic_distribution_Scaled.pdf")
savefig(p, "N1_N6_magic_distribution_Scaled.png")

# Mean, Median, Variance, Skewness, Kurtosis vs N

MostP       = Vector{Float64}()
Mean        = Vector{Float64}()
Variance    = Vector{Float64}()
Median      = Vector{Float64}()
Skewness    = Vector{Float64}()
Kurtosis    = Vector{Float64}()

for i in 1:6
    push!(MostP, MostProbable[i] ./ log2((2^i + 1)/2))
    push!(Mean, mean(Magic[i] ./ log2((2^i + 1)/2)) )
    push!(Variance, var(Magic[i] ./ log2((2^i + 1)/2)))
    push!(Median, median(Magic[i] ./ log2((2^i + 1)/2)))
    push!(Skewness, skewness(Magic[i] ./ log2((2^i + 1)/2)))
    push!(Kurtosis, kurtosis(Magic[i] ./ log2((2^i + 1)/2)))
end

p = plot(dpi=400)
x = 1:6
scatter!(x, Mean, label = "Mean")
scatter!(x, Median, label = "Median")
scatter!(x, MostP, label = L"max(P($M_2$))")
title!("Mean and Median for N = 1-6")
xlabel!(L"$N$",   labelfontsize = 20)
ylabel!(L"$\tilde{M}_2$",   labelfontsize = 20)
plot!(ylims=[0.55, 0.85])
plot!(framestyle=:box)
plot!(legendfontsize=10)

savefig(p, "N1_N6_mean_median_magic.pdf")
savefig(p, "N1_N6_mean_median_magic.png")

p = plot(dpi=400)
x = 1:6
scatter!(x, Variance, label = "Variance")
title!("Variance for N = 1-6")
xlabel!(L"$N$",   labelfontsize = 20)
#ylabel!(L"$\tilde{M}_2$",   labelfontsize = 20)
#plot!(ylims=[0.5, 1])
plot!(framestyle=:box)
plot!(legendfontsize=10)
plot!(yaxis=:log)
savefig(p, "N1_N6_variance_magic.pdf")
savefig(p, "N1_N6_variance_magic.png")

p = plot(dpi=400)
x = 1:6
scatter!(x, Skewness, label = "Skewness")
title!("Skewness for N = 1-6")
xlabel!(L"$N$",   labelfontsize = 20)
#ylabel!(L"$\tilde{M}_2$",   labelfontsize = 20)
#plot!(ylims=[0.5, 1])
plot!(framestyle=:box)
plot!(legendfontsize=10)

savefig(p, "N1_N6_skewness_magic.pdf")
savefig(p, "N1_N6_skewness_magic.png")

p = plot(dpi=400)
x = 1:6
scatter!(x, Kurtosis, label = "Skewness")
title!("Kurtosis for N = 1-6")
xlabel!(L"$N$",   labelfontsize = 20)
ylabel!(L"$\tilde{M}_2$",   labelfontsize = 20)
#plot!(ylims=[0.5, 1])
plot!(framestyle=:box)
plot!(legendfontsize=10)

savefig(p, "N1_N6_kurtosis_magic.pdf")
savefig(p, "N1_N6_kurtosis_magic.png")

# Entanglement and Magic
using KernelDensity
using StatsPlots

# N = 2
p = plot(dpi=400)
xlabel!(L"$\tilde{M}_2$", labelfontsize = 20)
ylabel!(L"$\tilde{S}$")
histogram2d!(Magic[2] ./ log2(5/2), map(first, data[2]["Entanglement"]), bins=(range(0, 1, 500), range(0.45, 0.55, 2)), show_empty_bins=false, normalise=:pdf)
k = kde((Magic[2] ./ log2(5/2), map(first, data[2]["Entanglement"])))
contour!(k, lw=1, levels=0:0.5:40)
plot!(xlims=[0, 1])
plot!(ylims=[0, 1])
title!(L"$\varrho(M_2, S)$ for N = 2")
plot!(framestyle=:box)
plot!(grid=false)
plot!(legendfontsize=10)

savefig(p, "N2_ms_dist.pdf")
savefig(p, "N2_ms_dist.png")

# Cuts along maximum magic and entropy:
h = fit(Histogram, (Magic[2] ./ log2(5/2), map(first, data[2]["Entanglement"])), (range(0,1,500), range(0, 1, 500)))

SCut = h.weights[:, 240:260]
MCut = h.weights[340:360, :]
p = plot(range(0, 1, 499), mean(SCut, dims=2) ./ sum(mean(SCut, dims=2)), lw=2, label="Max Entanglement cut")
plot!(range(0, 1, 499), mean(MCut, dims=1)[:] ./ sum(mean(MCut, dims=1)[:]), lw=2, label="Max Magic cut")
xlabel!(L"$\tilde{S}$ or $\tilde{M}$")
ylabel!(L"$\varrho(\tilde{S})$ at max($M$) or $\varrho(\tilde{M})$ at max($S$)")
savefig(p, "N2_sm_cuts.pdf")
savefig(p, "N2_sm_cuts.png")

# N = 3
N = 3
p = plot(dpi=400)

S = Vector{Float64}()
for i in 1:(2^20)
    push!(S, round(maximum(data[3]["Entanglement"][i][end]), digits=14))
end

histogram2d!(Magic[N] ./ log2((2^N + 1)/2), S, bins=(range(0, 1, 500), range(0, 1, 500)), show_empty_bins=false, normalise=:pdf)
k = kde((Magic[N] ./ log2((2^N + 1)/2), S))
maxIndex = findmax(k.density)
k.x[190]
contour!(k,lw=1, levels=100)

plot!(xlims=[0, 1])
plot!(ylims=[0, 1])

title!(L"$\varrho(M_2, S)$ for N = 3")
xlabel!(L"$\tilde{M}_2$", labelfontsize = 20)
ylabel!(L"$\tilde{S}$")

plot!(framestyle=:box)
plot!(grid=false)
plot!(legendfontsize=10)

savefig(p, "N3_ms_dist.pdf")
savefig(p, "N3_ms_dist.png")

# Cuts along maximum magic and entropy:
h = fit(Histogram, (Magic[N] ./ log2((2^N + 1)/2), S), (range(0,1,500), range(0, 1, 500)))

SCut = h.weights[:, (Int(0.72*500) - 10):(Int(0.72*500) + 10)]
MCut = h.weights[(Int(0.72*500) - 10):(Int(0.72*500) + 10), :]
p = plot(range(0, 1, 499), mean(SCut, dims=2) ./ sum(mean(SCut, dims=2)), lw=2, label="Max Entanglement cut")
plot!(range(0, 1, 499), mean(MCut, dims=1)[:] ./ sum(mean(MCut, dims=1)[:]), lw=2, label="Max Magic cut")
xlabel!(L"$\tilde{S}$ or $\tilde{M}$")
ylabel!(L"$\varrho(\tilde{S})$ at max($M$) or $\varrho(\tilde{M})$ at max($S$)")
savefig(p, "N$(N)_sm_cuts.pdf")
savefig(p, "N$(N)_sm_cuts.png")

# N = 4
N = 4
p = plot(dpi=400)

S = Vector{Float64}()
for i in 1:(2^20)
    push!(S, round(maximum(data[N]["Entanglement"][i][end]), digits=14))
end

S /= log2(4)
histogram2d!(Magic[N] ./ log2((2^N + 1)/2), S, bins=(range(0, 1, 500), range(0, 1, 500)), show_empty_bins=false, normalise=:pdf)
k = kde((Magic[N] ./ log2((2^N + 1)/2), S))
maxIndex = findmax(k.density)
k.x[198]
k.y[163]
contour!(k,lw=1, levels=10)

plot!(xlims=[0, 1])
plot!(ylims=[0, 1])

title!(L"$\varrho(M_2, S)$ for N = 4")
xlabel!(L"$\tilde{M}_2$", labelfontsize = 20)
ylabel!(L"$\tilde{S}$")

plot!(framestyle=:box)
plot!(grid=false)
plot!(legendfontsize=10)

savefig(p, "N4_ms_dist.pdf")
savefig(p, "N4_ms_dist.png")

# Cuts along maximum magic and entropy:
h = fit(Histogram, (Magic[N] ./ log2((2^N + 1)/2), S), (range(0,1,500), range(0, 1, 500)))

SCut = h.weights[:, (Int(0.68*500) - 10):(Int(0.68*500) + 10)]
MCut = h.weights[(Int(0.75*500) - 10):(Int(0.75*500) + 10), :]
p = plot(range(0, 1, 499), mean(SCut, dims=2) ./ sum(mean(SCut, dims=2)), lw=2, label="Max Entanglement cut")
plot!(range(0, 1, 499), mean(MCut, dims=1)[:] ./ sum(mean(MCut, dims=1)[:]), lw=2, label="Max Magic cut")
xlabel!(L"$\tilde{S}$ or $\tilde{M}$")
ylabel!(L"$\varrho(\tilde{S})$ at max($M$) or $\varrho(\tilde{M})$ at max($S$)")
savefig(p, "N$(N)_sm_cuts.pdf")
savefig(p, "N$(N)_sm_cuts.png")

# N = 5
N = 5
p = plot(dpi=400)

S = Vector{Float64}()

for i in 1:(2^20)
    push!(S, round((data[N]["Entanglement"][i][end]), digits=14))
end

S /= log2(4)
histogram2d!(Magic[N] ./ log2((2^N + 1)/2), S, bins=(range(0, 1, 500), range(0, 1, 500)), show_empty_bins=false, normalise=:pdf)
#
k = kde((Magic[N] ./ log2((2^N + 1)/2), S))
maxIndex = findmax(k.density)
k.x[194]
k.y[186]
#contour!(k,lw=1, levels=10)


plot!(xlims=[0, 1])
plot!(ylims=[0, 1])

title!(L"$\varrho(M_2, S)$ for N = 5")
xlabel!(L"$\tilde{M}_2$", labelfontsize = 20)
ylabel!(L"$\tilde{S}$")

plot!(framestyle=:box)
plot!(grid=false)
plot!(legendfontsize=10)

savefig(p, "N5_ms_dist.pdf")
savefig(p, "N5_ms_dist.png")

# Cuts along maximum magic and entropy:
h = fit(Histogram, (Magic[N] ./ log2((2^N + 1)/2), S), (range(0,1,500), range(0, 1, 500)))

SCut = h.weights[:, (Int(0.85*500) - 10):(Int(0.85*500) + 10)]
MCut = h.weights[(Int(0.78*500) - 10):(Int(0.78*500) + 10), :]
p = plot(range(0, 1, 499), mean(SCut, dims=2) ./ sum(mean(SCut, dims=2)), lw=2, label="Max Entanglement cut")
plot!(range(0, 1, 499), mean(MCut, dims=1)[:] ./ sum(mean(MCut, dims=1)[:]), lw=2, label="Max Magic cut")
xlabel!(L"$\tilde{S}$ or $\tilde{M}$")
ylabel!(L"$\varrho(\tilde{S})$ at max($M$) or $\varrho(\tilde{M})$ at max($S$)")
savefig(p, "N$(N)_sm_cuts.pdf")
savefig(p, "N$(N)_sm_cuts.png")


# N = 6
N = 6
p = plot(dpi=400)

S = data[6]["Svn_middle"]

S /= log2(8)
histogram2d!(Magic[N] ./ log2((2^N + 1)/2), S, bins=(range(0, 1, 500), range(0, 1, 500)), show_empty_bins=false, normalise=:pdf)
k = kde((Magic[N] ./ log2((2^N + 1)/2), S))
maxIndex = findmax(k.density)
k.x[191]
k.y[146]
contour!(k,lw=1, levels=10)


fit(Histogram, data[NIndex]["Magic"], bin_edges)
plot!(xlims=[0, 1])
plot!(ylims=[0, 1])

title!(L"$\varrho(M_2, S)$ for N = 6")
xlabel!(L"$\tilde{M}_2$", labelfontsize = 20)
ylabel!(L"$\tilde{S}$")

plot!(framestyle=:box)
plot!(grid=false)
plot!(legendfontsize=10)

savefig(p, "N6_ms_dist.pdf")
savefig(p, "N6_ms_dist.png")

# Cuts along maximum magic and entropy:
h = fit(Histogram, (Magic[N] ./ log2((2^N + 1)/2), S), (range(0,1,500), range(0, 1, 500)))

SCut = h.weights[:, (Int(0.85*500) - 10):(Int(0.85*500) + 10)]
MCut = h.weights[(Int(0.78*500) - 10):(Int(0.78*500) + 10), :]
p = plot(range(0, 1, 499), mean(SCut, dims=2) ./ sum(mean(SCut, dims=2)), lw=2, label="Max Entanglement cut")
plot!(range(0, 1, 499), mean(MCut, dims=1)[:] ./ sum(mean(MCut, dims=1)[:]), lw=2, label="Max Magic cut")
xlabel!(L"$\tilde{S}$ or $\tilde{M}$")
ylabel!(L"$\varrho(\tilde{S})$ at max($M$) or $\varrho(\tilde{M})$ at max($S$)")
savefig(p, "N$(N)_sm_cuts.pdf")
savefig(p, "N$(N)_sm_cuts.png")

# Covariance between entanglement in the middle and magic
Cov = Vector{Float64}()

for n in 2:6
    if n < 6
        S = Vector{Float64}()
        for i in 1:(2^20)
            push!(S, round(mean(data[n]["Entanglement"][i][:]), digits=14))
        end
        S /= log2(2^(floor(n/2)))
        push!(Cov, cov(Magic[n] ./ log((2^n + 1)/2), S))
    else
        push!(Cov, cov(Magic[n] ./ log((2^n + 1)/2), Svn_mean ./ log2(8)))
    end
end

p = plot(dpi=400)

scatter!(2:6, (Cov), label="Cov(N)")
plot!(yaxis=:log)
title!(L"$\langle M S \rangle - \langle M \rangle \langle S \rangle$")
xlabel!(L"$N$", labelfontsize = 20)
ylabel!(L"Covariance")

plot!(framestyle=:box)
plot!(grid=false)
plot!(legendfontsize=10)

savefig(p, "Cov_N_1_6_mean.pdf")
savefig(p, "Cov_N_1_6_mean.png")

# Covariance sample dependence for N = 5
Cov_Ns = Vector{Float64}()
Variance = Vector{Float64}()

S = Vector{Float64}()
for i in 1:(2^20)
    push!(S, round((data[5]["Entanglement"][i][end]), digits=14))
end

for i in 2:20
    println(i)
    push!(Cov_Ns, cov(Magic[5][1:2^i] ./ log2(33/2), S[1:2^i] ./ log2(4)))
end

p = plot(dpi=400)
scatter!(2:20, (Cov_Ns))

# Jackknife for the covariance N = 6
m = Magic[2]
s = data[2]["Svn_middle"]
sm = s .* m

JKCOV = Vector{Float64}()
for datalength in ProgressBar(5:10)
    M_means = Vector{Float64}()
    S_means = Vector{Float64}()
    SM_means = Vector{Float64}()
    for i in 1:2^
        push!(M_means, mean(m[[1:i-1; i+1:end]]))
        push!(S_means, mean(s[[1:i-1; i+1:end]]))
        push!(SM_means, mean(sm[[1:i-1; i+1:end]]))
    end
    push!(JKCOV, mean(SM_means) - mean(M_means)*mean(S_means))
end
covs = Vector{Float64}()
stds = Vector{Float64}()
for i in ProgressBar(1:20)
    push!(covs, cov(m[1:2^i], s[1:2^i]))
    push!(stds, std(m[1:2^i] .* s[1:2^i]))
end
scatter(1:20, abs.(covs), yaxis=:log)
scatter!(1:20, stds)
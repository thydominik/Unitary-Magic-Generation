using JLD2
using Statistics
using Plots
using DelimitedFiles
using ProgressBars
using StatsBase
# Figure 1 - Magic Distribution

PathN = Dict()
PathN[1] = "D:\\Data\\Random_Unitary_Magic_Generation\\N1\\Regular\\RegularUnitaryCircuitMagicSampled_N_1_Samples_1048576_Seed_1.jld2"
PathN[2] = "D:\\Data\\Random_Unitary_Magic_Generation\\N2\\Regular\\RegularUnitaryCircuitMagicSampled_N_2_Samples_1048576_Seed_1.jld2"
PathN[3] = "D:\\Data\\Random_Unitary_Magic_Generation\\N3\\Regular\\RegularUnitaryCircuitMagicSampled_N_3_Samples_1048576_Seed_1.jld2"
PathN[4] = "D:\\Data\\Random_Unitary_Magic_Generation\\N4\\Regular\\RegularUnitaryCircuitMagicSampled_N_4_Samples_1048576_Seed_1.jld2"
PathN[5] = "D:\\Data\\Random_Unitary_Magic_Generation\\N5\\Regular\\RegularUnitaryCircuitMagicSampled_N_5_Samples_1048576_Seed_1.jld2"
PathN[6] = "D:\\Data\\Random_Unitary_Magic_Generation\\N6\\Regular\\RegularUnitaryCircuitMagicSampled_N_6_Samples_1048576_MultiSeed_w_Ent.jld2"
PathN[7] = "D:\\Data\\Random_Unitary_Magic_Generation\\N7\\Regular\\RegularUnitaryCircuitMagicSampled_N_7_Samples_442368_MultiSeed_w_Ent.jld2"
PathN[8] = "D:\\Data\\Random_Unitary_Magic_Generation\\N8\\RegularUnitaryCircuitMagicSampled_N_8_Samples_425984.jld2"
PathN[9] = "D:\\Data\\Random_Unitary_Magic_Generation\\N9\\RegularUnitaryCircuitMagicSampled_N_9_Samples_49152.jld2"
PathN[10] = "D:\\Data\\Random_Unitary_Magic_Generation\\N10\\RegularUnitaryCircuitMagicSampled_N_10_Samples_4096_Seed_1.jld2"


Data = Dict()
for i in ProgressBar(8:10)
    Data[i] = JLD2.load(PathN[i])
end

#The data is in log() not log2()
Magic = Dict()
NMagic = Dict()

for i in 8:10
    if i < 8
        Magic[i] = -log2.(exp.(-(Data[i]["Magic"] .+ log(2^i)))) .- log2(2^i)
    else
        Magic[i] = Data[i]["Magic"] #-log2.(exp.(-(Data[i]["Magic"] .+ log(2^i)))) .- log2(2^i)
    end
    NMagic[i] = Magic[i] ./ log2((2^i + 1)/2)
    #Magic[i] = -log2.(-exp.((Data[i]["Magic"] .+ log(2^i)))) .- log2(2^i)
end 


Mean        = Vector()
Meanstd     = Vector()
Variance    = Vector()
Variancestd = Vector()
Median      = Vector()
Medianstd   = Vector()
Skewness    = Vector()
Skewnessstd = Vector()
Kurtosis    = Vector()
Kurtosisstd = Vector()

for i in 1:10
    data = NMagic[i];
    push!(Mean, mean(data))
    push!(Variance, var(data))
    push!(Median, median(data))
    push!(Skewness, skewness(data))
    push!(Kurtosis, kurtosis(data))
    
    results = bootstrap_statistics(data, 5000)

    push!(Meanstd, std(results[1]))
    push!(Variancestd, std(results[2]))
    push!(Medianstd, std(results[3]))
    push!(Skewnessstd, std(results[4]))
    push!(Kurtosisstd, std(results[5]))
end

N = 1:8
writedlm("statistical_Variables.txt", hcat(N, Mean, Meanstd, Median, Medianstd, Variance, Variancestd, Skewness, Skewnessstd, Kurtosis, Kurtosisstd))


# get middle entanglement
SVN = Dict()
for i in 8:10
    println(i)
    if i < 10
        SVN[i] = map(last, Data[i]["Svn"])
    else
        SVN[i] = map(last, Data[i]["Entanglement"])
    end
end

SVN[6] = Data[6]["Svn_middle"]
SVN[7] = Data[7]["Svn_middle"]
SVN[8] = Data[8]["Svn_max"]

COV = Vector()
COVstd = Vector()
for i in 8:10
    M = NMagic[i];
    SE = SVN[i] ./ log(2^floor(i/2));
    result = covariance_with_error(M, SE, n_bootstrap=1000)
    push!(COV, cov(NMagic[i], SE))
    push!(COVstd, abs(result[2]))

end

writedlm("Covariance.txt", hcat(2:8, abs.(COV), abs.(COVstd)./2))

function covariance_with_error(data1::Vector{Float64}, data2::Vector{Float64}; n_bootstrap::Int=1000)
    # Calculate the covariance
    cov_value = cov(data1, data2)
    
    # Bootstrap function for covariance
    function bootstrap_cov()
        indices = rand(1:length(data1), length(data1))
        return cov(data1[indices], data2[indices])
    end
    
    # Perform bootstrap
    bootstrap_covs = [bootstrap_cov() for _ in ProgressBar(1:n_bootstrap)]
    
    # Calculate the standard error
    cov_error = std(bootstrap_covs)
    
    return cov_value, cov_error
end



# S and M distribution
h2 = Dict()
for i in 8:10
    SE = SVN[i] ./ log2(2^floor(i/2));
    h2[i] = fit(Histogram, (NMagic[i], SE), (range(0, 1, 500), range(0, 1, 500)))
    h2[i] = StatsBase.normalize(h2[i], mode=:pdf)

    writedlm("N$(i)_S_and_M_dist.txt", h2[i].weights)
end
writedlm("s_and_m_distribiution_x_and_y.txt", hcat(h2[2].edges[1], h2[2].edges[2]))

plot(h2[10])














using StatsBase

# Assuming you have your data in an array called 'data'

data = NMagic[7]  # Example: 1000 random numbers from standard normal distribution

histogram(data, nbins=100)
# Calculate basic statistics
mean_value = mean(data)
variance_value = var(data)
median_value = median(data)
kurtosis_value = kurtosis(data)
skewness_value = skewness(data)

using Statistics

function bootstrap_statistics(data, n_bootstrap=10000)
    n = length(data)
    bootstrap_means = zeros(n_bootstrap)
    bootstrap_variances = zeros(n_bootstrap)
    bootstrap_medians = zeros(n_bootstrap)
    bootstrap_kurtoses = zeros(n_bootstrap)
    bootstrap_skewnesses = zeros(n_bootstrap)
    
    for i in ProgressBar(1:n_bootstrap)
        sample = data[rand(1:n, n)]
        bootstrap_means[i] = mean(sample)
        bootstrap_variances[i] = var(sample)
        bootstrap_medians[i] = median(sample)
        bootstrap_kurtoses[i] = kurtosis(sample)
        bootstrap_skewnesses[i] = skewness(sample)
    end
    
    return bootstrap_means, bootstrap_variances, bootstrap_medians, bootstrap_kurtoses, bootstrap_skewnesses
end

# Perform bootstrap
bootstrap_results = bootstrap_statistics(data)

# Calculate standard errors
se_mean = std(bootstrap_results[1])
se_variance = std(bootstrap_results[2])
se_median = std(bootstrap_results[3])
se_kurtosis = std(bootstrap_results[4])
se_skewness = std(bootstrap_results[5])

println("Mean: $(mean_value) ± $(se_mean)")
println("Variance: $(variance_value) ± $(se_variance)")
println("Median: $(median_value) ± $(se_median)")
println("Kurtosis: $(kurtosis_value) ± $(se_kurtosis)")
println("Skewness: $(skewness_value) ± $(se_skewness)")

using Plots

statistics = [mean_value, variance_value, median_value, kurtosis_value, skewness_value]
errors = [se_mean, se_variance, se_median, se_kurtosis, se_skewness]
labels = ["Mean", "Variance", "Median", "Kurtosis", "Skewness"]

plot(1:5, statistics, yerr=errors, label="", xlabel="Statistic", ylabel="Value",
     xticks=(1:5, labels), marker=:circle)


mean(NMagic[1])
mean(h[1].weights)

h = Dict()
bins = Dict()
for i in 1:8
    bins[i] = range(0, log2((2^i + 1)/2), 1000)
    h[i] = fit(Histogram, Data[i]["Magic"], bins[i])
    h[i] = StatsBase.normalize(h[i], mode=:pdf)
end

p = plot(dpi=300)
for i in 1:8
    plot!(collect(bins[i])[1:end-1], h[i].weights, lw=2)
end
plot!()
savefig(p, "N1_8_histogram.png")
for i in 1:8
    savedata = hcat(collect(h[i].edges[1])[1:end-1], h[i].weights)
    writedlm("N$(i)_histogram.txt", savedata)
end

Nh = Dict()
bins = Dict()
p=plot(dpi=300)
for i in 1:8
    bins[i] = range(0, 1, 1000)
    Nh[i] = fit(Histogram, NMagic[i], bins[i])
    Nh[i] = StatsBase.normalize(Nh[i], mode=:pdf)
    plot!(collect(bins[i])[1:end-1], Nh[i].weights, lw=2)
end
plot!()
savefig(p, "N1_8_norm_histogram.png")

for i in 1:8
    savedata = hcat(collect(bins[i])[1:end-1], Nh[i].weights)
    writedlm("N$(i)_norm_histogram.txt", savedata)
end

# Depth dependence N = 5

Path = "D:\\Data\\Random_Unitary_Magic_Generation\\N5\\BrickWall\\"
Data = Dict()
for i in 1:6
    Data[i] = JLD2.load(Path * "BWUnitaryCircuitMagicSampled_N_5_D_$(i)_Samples_1048576_Seed_1.jld2")
end

#The data is in log() not log2()
Magic = Dict()
NMagic = Dict()
for i in 1:6
    Magic[i] = -log2.(exp.(-(Data[i]["Magic"] .+ log(2^5)))) .- log2(2^5)
    NMagic[i] = Magic[i] ./ log2((2^5 + 1)/2)
    #Magic[i] = -log2.(-exp.((Data[i]["Magic"] .+ log(2^i)))) .- log2(2^i)
end

h = Dict()
bins = Dict()
for i in 1:6
    bins[i] = range(0, 1, 1000)
    h[i] = fit(Histogram, NMagic[i], bins[i])
    h[i] = StatsBase.normalize(h[i], mode=:pdf)
end
p = plot(dpi=300)
for i in 1:6
    plot!(collect(bins[i])[1:end-1], h[i].weights, lw=2)
end
plot!()

for i in 1:6
    savedata = hcat(collect(bins[i])[1:end-1], h[i].weights)
    writedlm("N5_D$(i)_norm_histogram.txt", savedata)
end







Data = Dict()
for i in [7, 8, 9, 10, 15]
    Data[i] = JLD2.load("BWUnitaryCircuitMagicSampled_N_5_D_$(i)_Samples_1048576_Seed_1.jld2")
end

#The data is in log() not log2()
Magic = Dict()
NMagic = Dict()
for i in [7, 8, 9, 10, 15]
    Magic[i]    = Data[i]["Magic"]
    #Magic[i]    = -log2.(exp.(-(Data[i]["Magic"] .+ log(2^5)))) .- log2(2^5)
    NMagic[i]   = Magic[i] ./ log2((2^5 + 1)/2)
    
end

h = Dict()
bins = Dict()
for i in [7, 8, 9, 10, 15]
    bins[i] = range(0, 1, 1000)
    h[i] = fit(Histogram, NMagic[i], bins[i])
    h[i] = StatsBase.normalize(h[i], mode=:pdf)
end
p = plot(dpi=300)
for i in [7, 8, 9, 10, 15]
    plot!(collect(bins[i])[1:end-1], h[i].weights, lw=2)
end
plot!()

for i in [7, 8, 9, 10, 15]
    savedata = hcat(collect(bins[i])[1:end-1], h[i].weights)
    writedlm("N5_D$(i)_norm_histogram.txt", savedata)
end







using DelimitedFiles

Data = Dict()
for i in 1:8
    Data[i] = readdlm("N$(i)_histogram.txt")
end

p = plot(dpi=300)
for i in 1:8
    println(log2((2^i + 1)/2) / log((2^i +1)/2))
    plot!(Data[i][:, 1] .* log2((2^i + 1)/2) / log((2^i +1)/2), Data[i][:, 2])
    writedlm("New_N$(i)_histogram.txt", hcat(Data[i][:, 1].* log2((2^i + 1)/2) / log((2^i +1)/2), Data[i][:, 2] ./ 1.4426950408889634))
end
plot!(xlims=(0, log2((2^8 + 1)/2)))

using DelimitedFiles

Data = Dict()
for i in 1:8
    Data[i] = readdlm("New_N$(i)_histogram.txt")
end

p = plot(dpi=300)
for i in 1:8
    plot!(Data[i][:, 1], Data[i][:, 2]  )
end
plot!(xlims=(0, log2((2^8 + 1)/2)))

Integr = 0
for i in 2:length(Data[5][:, 2])
    Integr += (Data[5][i, 1] - Data[5][i - 1, 1]) * 0.5 * (Data[5][i, 2] + Data[5][i - 1, 2])
end
println(Integr)
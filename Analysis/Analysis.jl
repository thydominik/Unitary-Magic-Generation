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




using JLD2
using Statistics
using Plots
using DelimitedFiles
using ProgressBars
using StatsBase
# Figure 1 - Magic Distribution

PathN = Dict()
PathN[1]    = "D:\\Data\\Random_Unitary_Magic_Generation\\Simplify\\RegularUnitaryCircuitMagicSampled_N_1_Samples_1048576.jld2"
PathN[2]    = "D:\\Data\\Random_Unitary_Magic_Generation\\Simplify\\RegularUnitaryCircuitMagicSampled_N_2_Samples_1048576.jld2"
PathN[3]    = "D:\\Data\\Random_Unitary_Magic_Generation\\Simplify\\RegularUnitaryCircuitMagicSampled_N_3_Samples_1048576.jld2"
PathN[4]    = "D:\\Data\\Random_Unitary_Magic_Generation\\Simplify\\RegularUnitaryCircuitMagicSampled_N_4_Samples_1048576.jld2"
PathN[5]    = "D:\\Data\\Random_Unitary_Magic_Generation\\Simplify\\RegularUnitaryCircuitMagicSampled_N_5_Samples_1048576.jld2"
PathN[6]    = "D:\\Data\\Random_Unitary_Magic_Generation\\Simplify\\RegularUnitaryCircuitMagicSampled_N_6_Samples_1048576.jld2"
PathN[7]    = "D:\\Data\\Random_Unitary_Magic_Generation\\Simplify\\RegularUnitaryCircuitMagicSampled_N_7_Samples_442368.jld2"
PathN[8]    = "D:\\Data\\Random_Unitary_Magic_Generation\\Simplify\\RegularUnitaryCircuitMagicSampled_N_8_Samples_1048576.jld2"
PathN[9]    = "D:\\Data\\Random_Unitary_Magic_Generation\\Simplify\\RegularUnitaryCircuitMagicSampled_N_9_Samples_196608.jld2"
PathN[10]   = "D:\\Data\\Random_Unitary_Magic_Generation\\Simplify\\RegularUnitaryCircuitMagicSampled_N_10_Samples_20480.jld2"


Data = Dict()
for i in ProgressBar(1:10)
    Data[i] = JLD2.load(PathN[i])
end

#The data is in log() not log2()
Magic = Dict()
NMagic = Dict()

for i in 1:10
    if i < 8
        Magic[i] = -log2.(exp.(-(Data[i]["Magic"] .+ log(2^i)))) .- log2(2^i)
    else
        Magic[i] = Data[i]["Magic"] 
    end
    NMagic[i] = Magic[i] ./ log2((2^i + 1)/2)
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
    data = Magic[i];
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

N = 1:10
writedlm("statistical_Variables.txt", hcat(N, Mean, Meanstd, Median, Medianstd, Variance, Variancestd, Skewness, Skewnessstd, Kurtosis, Kurtosisstd))


# get middle entanglement
SVN = Dict()
for i in 2:10
    SVN[i] = Data[i]["Svn"]
end

# Covariance calculation
begin
    COV = Vector()
    for i in 2:9
        M   = Magic[i];
        SE  = SVN[i] # ./ log2(2^floor(i/2));
        push!(COV, cov(M, SE))
    end

println(COV)
COV[6] /= 2.37
COV[8] /= 5.333
writedlm("Covariance_M_S.txt", hcat(2:9, abs.(COV)))
end

# ----------------------------------------------
# Covariance error calculation

COVstd  = Vector()
for i in 2:10
    M   = Magic[i];
    SE  = SVN[i] # ./ log2(2^floor(i/2));
    result = covariance_with_error(M, SE, n_bootstrap=5000)
    push!(COV, cov(M, SE))
    push!(COVstd, abs(result[2]))
end
COVstd[6] /= 2.37
COVstd[7] /= 2
COVstd[8] *= 2
COVstd /= 2
writedlm("Covariance_error.txt", hcat(2:10, abs.(COVstd)))

# ----------------------------------------------

# S and M distribution
h2 = Dict()
for i in 2:10
    SE = SVN[i] #./ log2(2^floor(i/2));
    h2[i] = fit(Histogram, (Magic[i], SE), (range(minimum(Magic[i]), maximum(Magic[i]), 500), range(minimum(SE), maximum(SE), 500)))
    h2[i] = StatsBase.normalize(h2[i], mode=:pdf)

    writedlm("N$(i)_S_and_M_dist.txt", h2[i].weights)
    writedlm("N$(i)_S_and_M_dist_bins.txt", hcat(h2[i].edges[1], h2[i].edges[2]))
end

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




# Entropy histograms

h = Dict()
bins = Dict()
for i in 2:10
    binscale = Int(round(log2(2^floor(i/2)) / 1) * 250)
    bins[i] = range(0, log2(2^floor(i/2)), binscale)
    h[i] = fit(Histogram, SVN[i], bins[i])
    h[i] = StatsBase.normalize(h[i], mode=:pdf)
end

EntropyMean        = Vector()
EntropyMeanstd     = Vector()
EntropyVariance    = Vector()
EntropyVariancestd = Vector()
EntropyMedian      = Vector()
EntropyMedianstd   = Vector()
EntropySkewness    = Vector()
EntropySkewnessstd = Vector()
EntropyKurtosis    = Vector()
EntropyKurtosisstd = Vector()

for i in 2:10
    data = SVN[i];
    push!(EntropyMean, mean(data))
    push!(EntropyVariance, var(data))
    push!(EntropyMedian, median(data))
    push!(EntropySkewness, skewness(data))
    push!(EntropyKurtosis, kurtosis(data))
    
    results = bootstrap_statistics(data, 5000)

    push!(EntropyMeanstd, std(results[1]))
    push!(EntropyVariancestd, std(results[2]))
    push!(EntropyMedianstd, std(results[3]))
    push!(EntropySkewnessstd, std(results[4]))
    push!(EntropyKurtosisstd, std(results[5]))
end

writedlm("Entropy_Statistical_Variables.txt", hcat(2:10, EntropyMean, EntropyMeanstd, EntropyMedian, EntropyMedianstd, EntropyVariance, EntropyVariancestd, EntropySkewness, EntropySkewnessstd, EntropyKurtosis, EntropyKurtosisstd))


p = plot(dpi=300)
for i in 2:10
    plot!(collect(bins[i])[1:end-1], h[i].weights, lw=2)
end
plot!()
savefig(p, "N1_8_histogram.png")
for i in 2:10
    savedata = hcat(collect(h[i].edges[1])[1:end-1], h[i].weights)
    writedlm("N$(i)_entanglement_histogram.txt", savedata)
end

h = Dict()
bins = Dict()
for i in 1:10
    binscale = Int(round(log2((2^i + 1)/2) / log2(3/2)) * 500)
    bins[i] = range(0, log2((2^i + 1)/2), binscale)
    h[i] = fit(Histogram, Magic[i], bins[i])
    h[i] = StatsBase.normalize(h[i], mode=:pdf)
end


Nh = Dict()
bins = Dict()
p=plot(dpi=300)
for i in 1:8
    bins[i] = range(0, 1, 1000)
    Nh[i] = fit(Histogram, Magic[i], bins[i])
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
using JLD2
using Statistics
using Plots
using DelimitedFiles
# Figure 1 - Magic Distribution

PathN = Dict()
PathN[1] = "D:\\Data\\Random_Unitary_Magic_Generation\\N1\\Regular\\RegularUnitaryCircuitMagicSampled_N_1_Samples_1048576_Seed_1.jld2"
PathN[2] = "D:\\Data\\Random_Unitary_Magic_Generation\\N2\\Regular\\RegularUnitaryCircuitMagicSampled_N_2_Samples_1048576_Seed_1.jld2"
PathN[3] = "D:\\Data\\Random_Unitary_Magic_Generation\\N3\\Regular\\RegularUnitaryCircuitMagicSampled_N_3_Samples_1048576_Seed_1.jld2"
PathN[4] = "D:\\Data\\Random_Unitary_Magic_Generation\\N4\\Regular\\RegularUnitaryCircuitMagicSampled_N_4_Samples_1048576_Seed_1.jld2"
PathN[5] = "D:\\Data\\Random_Unitary_Magic_Generation\\N5\\Regular\\RegularUnitaryCircuitMagicSampled_N_5_Samples_1048576_Seed_1.jld2"
PathN[6] = "D:\\Data\\Random_Unitary_Magic_Generation\\N6\\Regular\\RegularUnitaryCircuitMagicSampled_N_6_Samples_1048576_MultiSeed_w_Ent.jld2"
PathN[7] = "D:\\Data\\Random_Unitary_Magic_Generation\\N7\\Regular\\RegularUnitaryCircuitMagicSampled_N_7_Samples_442368_MultiSeed_w_Ent.jld2"
PathN[8] = "D:\\Data\\Random_Unitary_Magic_Generation\\N8\\RegularUnitaryCircuitMagicSampled_N_8_Samples_81920_Seed_1_w_Ent.jld2"

Data = Dict()
for i in 1:8
    Data[i] = JLD2.load(PathN[i])
end

#The data is in log() not log2()
Magic = Dict()
NMagic = Dict()
for i in 1:8
    Magic[i] = -log2.(exp.(-(Data[i]["Magic"] .+ log(2^i)))) .- log2(2^i)
    NMagic[i] = Magic[i] ./ log2((2^i + 1)/2)
    #Magic[i] = -log2.(-exp.((Data[i]["Magic"] .+ log(2^i)))) .- log2(2^i)
end 

h = Dict()
bins = Dict()
for i in 1:8
    bins[i] = range(0, log((2^i + 1)/2), 1000)
    h[i] = fit(Histogram, Data[i]["Magic"], bins[i])
    h[i] = StatsBase.normalize(h[i], mode=:pdf)
end

plot(dpi=300)
for i in 1:8
    plot!(collect(bins[i])[1:end-1], h[i].weights, lw=2)
end

plot(dpi=300)
for i in 1:8
    histogram!(NMagic[i], bins=range(0,1, 2000), normalize=:pdf)
end
plot!()
histogram(NMagic[7])
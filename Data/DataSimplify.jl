using JLD2
using Plots
using StatsBase
using ProgressBars

NoSamples = 442368
NoSeeds = 1
N = 7

FilePath = "D:\\Data\\Random_Unitary_Magic_Generation\\N$(N)\\Regular\\"
Magic       = Vector{Float64}()
#Svn_max     = Vector{Float64}()
#Svn_middle  = Vector{Float64}()
#Svn_mean    = Vector{Float64}()
Svn         = Vector{Float64}()

D = []
for seedIndex in ProgressBar(1:NoSeeds)
    datafile = "RegularUnitaryCircuitMagicSampled_N_6_Samples_1048576_MultiSeed_w_Ent"
    
    #datafile = "BWUnitaryCircuitMagicSampled_N_7_D_$(Depth)_Samples_32768_Seed_$(seedIndex).jld2";

    #D = load(FilePath * datafile);
    D = load("D:\\Data\\Random_Unitary_Magic_Generation\\N7\\Regular\\RegularUnitaryCircuitMagicSampled_N_7_Samples_442368_MultiSeed_w_Ent.jld2")
    for sampleIndex in ProgressBar(1:NoSamples)
        #push!(Svn_max, maximum(D["Entanglement"][sampleIndex]))
        #push!(Svn_mean, mean(D["Entanglement"][sampleIndex]))
        #push!(Svn_middle, D["Entanglement"][sampleIndex][end])
        #push!(Svn, D["Entanglement"][sampleIndex][end])
        push!(Svn, D["Svn_middle"][sampleIndex])
        push!(Magic, D["Magic"][sampleIndex])
    end
end
    #saveFileName = "BWUnitaryCircuitMagicSampled_N_7_D_$(Depth)_Samples_$(10*32768)_MultiSeed_w_Ent.jld2"
saveFileName = "RegularUnitaryCircuitMagicSampled_N_$(N)_Samples_$(NoSeeds * NoSamples).jld2"
@save saveFileName Magic Svn #Svn_max Svn_mean Svn_middle

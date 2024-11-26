using JLD2
using Plots
using StatsBase
using ProgressBars
FilePath = "D:\\Data\\Random_Unitary_Magic_Generation\\N6\\Regular\\"

Magic       = Vector{Float64}()
Svn_max     = Vector{Float64}()
Svn_middle  = Vector{Float64}()
Svn_mean    = Vector{Float64}()

for seedIndex in ProgressBar(1:32)
    datafile = "RegularUnitaryCircuitMagicSampled_N_6_Samples_32768_Seed_$(seedIndex)_w_Ent.jld2"
    D = load(FilePath * datafile)
    for sampleIndex in 1:32768
        push!(Svn_max, maximum(D["Entanglement"][sampleIndex]))
        push!(Svn_mean, mean(D["Entanglement"][sampleIndex]))
        push!(Svn_middle, D["Entanglement"][sampleIndex][end])
        push!(Magic, D["Magic"][sampleIndex])
    end
end

saveFileName = "RegularUnitaryCircuitMagicSampled_N_6_Samples_1048576_MultiSeed_w_Ent.jld2"
@save saveFileName Magic Svn_max Svn_mean Svn_middle

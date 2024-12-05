using JLD2
using Plots
using StatsBase
using ProgressBars
FilePath = "D:\\Data\\Random_Unitary_Magic_Generation\\N8\\"
Magic       = Vector{Float64}()
Svn_max     = Vector{Float64}()
Svn_middle  = Vector{Float64}()
Svn_mean    = Vector{Float64}()


    for seedIndex in ProgressBar(1:9)

        datafile = "RegularUnitaryCircuitMagicSampled_N_8_Samples_8192_Seed_$(seedIndex)_w_Ent.jld2"
        #datafile = "BWUnitaryCircuitMagicSampled_N_7_D_$(Depth)_Samples_32768_Seed_$(seedIndex).jld2";
        D = load(FilePath * datafile);
        for sampleIndex in 1:8192
            push!(Svn_max, maximum(D["Entanglement"][sampleIndex]))
            push!(Svn_mean, mean(D["Entanglement"][sampleIndex]))
            push!(Svn_middle, D["Entanglement"][sampleIndex][end])
            push!(Magic, D["Magic"][sampleIndex])
        end
    end
    #saveFileName = "BWUnitaryCircuitMagicSampled_N_7_D_$(Depth)_Samples_$(10*32768)_MultiSeed_w_Ent.jld2"
    saveFileName = "RegularUnitaryCircuitMagicSampled_N_8_Samples_$(10*8192)_Seed_1_w_Ent.jld2"
    @save saveFileName Magic Svn_max Svn_mean Svn_middle


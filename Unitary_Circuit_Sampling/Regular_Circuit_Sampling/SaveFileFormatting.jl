using JLD2
using MAT

for N in 2:5
    fname1 = "RegularUnitaryCircuitMagicSampled_N_$(N)_Samples_1048576_Seed_1.jdl2"
    data = JLD2.load(fname1)
    fname2 = "RegularUnitaryCircuitMagicSampled_N_$(N)_Samples_1048576_Seed_1.mat"
    matwrite(fname2, data)
end

data = JLD2.load("RegularUnitaryCircuitMagicSampled_N_2_Samples_1048576_Seed_1.jdl2")
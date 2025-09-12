using Random
using ProgressBars
using JLD2
using Base.Threads
using MAT

current_dir = @__DIR__

# Unitary matrices
filepath = joinpath(current_dir, "..", "Modules", "Magic.jl")
include(filepath)
# Magic
filepath = joinpath(current_dir, "..", "Modules", "Random_Unitaries.jl")
include(filepath)
# Measure_Entanglement
filepath = joinpath(current_dir, "..", "Modules", "Entanglement.jl")
include(filepath)

using .Measure_Entanglement
using .Measure_Magic
using .Random_Unitary_Generation

Seed = 1
Random.seed!(Seed)
No_Samples = 2^21
No_Qubits = 2


# Generate all 2-qubit Clifford gates
Cliffords, labels = Generate_All_2_Qubit_Clifford_Gates()

# Start with the 2-qubit |00⟩ state
zero_state = [1.0 + 0im, 0.0 + 0im, 0.0 + 0im, 0.0 + 0im]

# Apply each Clifford to |00⟩ and collect resulting states
clifford_states = [C * zero_state for C in Cliffords]

# Find unique states (up to a global phase)
function unique_states(states)
    unique_vecs = []
    for v in states
        is_unique = true
        for u in unique_vecs
            # Check if vectors are equal up to a global phase
            phase = v[findfirst(!iszero, v)] / u[findfirst(!iszero, u)]
            if all(abs.(v .- phase * u) .< 1e-8)
                is_unique = false
                break
            end
        end
        if is_unique
            push!(unique_vecs, v)
        end
    end
    return unique_vecs
end

Psi = unique_states([round.(v; digits=10) for v in clifford_states])
No_Stab = length(Psi)






Psi1 = hcat(Psi...)'  # Convert to matrix with each state as a row
# Generate all 2-qubit stabilizer product states
No_Stab = No_Stab

Psi = Array{ComplexF64,2}(undef, No_Stab, 4)
idx = 1


Strings = Measure_Magic.GenerateAllPauliStrings(2)
PauliOperators = Measure_Magic.PauliOperatorList(Strings, 2)
Magic = zeros(No_Samples, No_Stab)

for i in ProgressBar(1:No_Samples)
    #M2 = Vector{Float64}()
    U = Random_Unitary_Generation.Generate_Regular_Unitary_Circuit(No_Qubits)
    for stabstate in 1:No_Stab
        Psi_0 = Psi1[stabstate, :]
        
        State = U * Psi_0
        #push!(M2, Measure_Magic.MeasureMagic_Pure(State, PauliOperators, 2)[1])
        Magic[i, stabstate] = Measure_Magic.MeasureMagic_Pure(State, PauliOperators, 2)[1]
    end
    #push!(Magic, sum(M2)/No_Stab)
end

fname = "NonStabilisingPower_of_Unitary_N_2_Samples_$(No_Samples)_Seed_$(Seed).mat"
matwrite(fname, Dict("Magic" => Magic))
fname = "UnitaryMagicSampling_AllStabs_N_2_Samples_$(No_Samples).jld2"
@save fname Magic

using Plots
using Statistics

plot(xlims=(0, log2(5/2)), ylims=(10^-3, 15))
histogram!(Magic[:, 3], bins=range(0, log2(5/2), 1000), normalize=:pdf)

NSP = zeros(No_Samples)
for i in 1:No_Samples
    NSP[i] = sum(Magic[i, :]) / No_Stab
end
histogram!(NSP, bins=range(0, log2(5/2), 1000), normalize=:pdf)
plot!(yaxis=:log10)

mean(Magic)
mean(Magic[:, 1])
mean(NSP) - mean(Magic[:, 2])




function Generate_All_2_Qubit_Clifford_Gates()
    # Identity
    Id = [1 0; 0 1]
    # Pauli Gates
    PX = [0 1; 1 0]
    PY = [0 -1im; 1im 0]
    PZ = [1 0; 0 -1]
    # Hadamard
    H = 1/sqrt(2) * [1 1; 1 -1]
    # Phase gate
    S = [1 0; 0 1im]
    
    # Two site Clifford qubit gates (11520):
    Clifford_Gates_2_Site = []
    
    # Hadamard "group" and its labels
    Hadamard_group = [Id, H]
    Hadamard_group_labels = ["I", "H"]
    
    # Phase shift "group" and its labels
    Phaseshift_group = [Id, H * S, S * H]
    Phaseshift_group_labels = ["I", "HS", "SH"]
    
    # "States" "group" and its labels
    Pauli_group = [Id, PX, PY, PZ]
    Pauli_group_labels = ["I", "X", "Y", "Z"]
    
    # Gates created from the 1 qubit clifford gates (2, 2, 3, 3, 4, 4)
    Clifford_Labels = String[]
    for hj1 in 1:length(Hadamard_group), hj2 in 1:length(Hadamard_group),
        vj1 in 1:length(Phaseshift_group), vj2 in 1:length(Phaseshift_group),
        pj1 in 1:length(Pauli_group), pj2 in 1:length(Pauli_group)
        
        push!(Clifford_Gates_2_Site, kron(Hadamard_group[hj1], Hadamard_group[hj2]) * 
                                     kron(Phaseshift_group[vj1], Phaseshift_group[vj2]) * 
                                     kron(Pauli_group[pj1], Pauli_group[pj2]))
        push!(Clifford_Labels, "($(Hadamard_group_labels[hj1])⊗$(Hadamard_group_labels[hj2]))*" *
                               "($(Phaseshift_group_labels[vj1])⊗$(Phaseshift_group_labels[vj2]))*" *
                               "($(Pauli_group_labels[pj1])⊗$(Pauli_group_labels[pj2]))")
    end
    
    # CNOT class contains 5184 elements = 24^2*3^2 or (2, 2, 3, 3, 3, 3, 4, 4)
    CNOT = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]
    
    for hj1 in 1:length(Hadamard_group), hj2 in 1:length(Hadamard_group),
        vj1 in 1:length(Phaseshift_group), vj2 in 1:length(Phaseshift_group),
        vj3 in 1:length(Phaseshift_group), vj4 in 1:length(Phaseshift_group),
        pj1 in 1:length(Pauli_group), pj2 in 1:length(Pauli_group)
        
        push!(Clifford_Gates_2_Site, kron(Hadamard_group[hj1], Hadamard_group[hj2]) * 
                                     kron(Phaseshift_group[vj1], Phaseshift_group[vj2]) * 
                                     CNOT * 
                                     kron(Phaseshift_group[vj3], Phaseshift_group[vj4]) * 
                                     kron(Pauli_group[pj1], Pauli_group[pj2]))
        push!(Clifford_Labels, "($(Hadamard_group_labels[hj1])⊗$(Hadamard_group_labels[hj2]))*" *
                               "($(Phaseshift_group_labels[vj1])⊗$(Phaseshift_group_labels[vj2]))*" *
                               "CNOT*" *
                               "($(Phaseshift_group_labels[vj3])⊗$(Phaseshift_group_labels[vj4]))*" *
                               "($(Pauli_group_labels[pj1])⊗$(Pauli_group_labels[pj2]))")
    end
    
    # class 3 elements with SWAP and CNOT class contains 5184 elements
    SWAP = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
    
    for hj1 in 1:length(Hadamard_group), hj2 in 1:length(Hadamard_group),
        vj1 in 1:length(Phaseshift_group), vj2 in 1:length(Phaseshift_group),
        vj3 in 1:length(Phaseshift_group), vj4 in 1:length(Phaseshift_group),
        pj1 in 1:length(Pauli_group), pj2 in 1:length(Pauli_group)
        
        push!(Clifford_Gates_2_Site, kron(Hadamard_group[hj1], Hadamard_group[hj2]) * 
                                     kron(Phaseshift_group[vj1], Phaseshift_group[vj2]) * 
                                     SWAP * CNOT * 
                                     kron(Phaseshift_group[vj3], Phaseshift_group[vj4]) * 
                                     kron(Pauli_group[pj1], Pauli_group[pj2]))
        push!(Clifford_Labels, "($(Hadamard_group_labels[hj1])⊗$(Hadamard_group_labels[hj2]))*" *
                               "($(Phaseshift_group_labels[vj1])⊗$(Phaseshift_group_labels[vj2]))*" *
                               "SWAP*CNOT*" *
                               "($(Phaseshift_group_labels[vj3])⊗$(Phaseshift_group_labels[vj4]))*" *
                               "($(Pauli_group_labels[pj1])⊗$(Pauli_group_labels[pj2]))")
    end
    
    # class 4 elements contains 576 elements (2, 2, 3, 3, 4, 4)
    for hj1 in 1:length(Hadamard_group), hj2 in 1:length(Hadamard_group),
        vj1 in 1:length(Phaseshift_group), vj2 in 1:length(Phaseshift_group),
        pj1 in 1:length(Pauli_group), pj2 in 1:length(Pauli_group)
        
        push!(Clifford_Gates_2_Site, kron(Hadamard_group[hj1], Hadamard_group[hj2]) * 
                                     kron(Phaseshift_group[vj1], Phaseshift_group[vj2]) * 
                                     SWAP * 
                                     kron(Pauli_group[pj1], Pauli_group[pj2]))
        push!(Clifford_Labels, "($(Hadamard_group_labels[hj1])⊗$(Hadamard_group_labels[hj2]))*" *
                               "($(Phaseshift_group_labels[vj1])⊗$(Phaseshift_group_labels[vj2]))*" *
                               "SWAP*" *
                               "($(Pauli_group_labels[pj1])⊗$(Pauli_group_labels[pj2]))")
    end
    
    return Clifford_Gates_2_Site, Clifford_Labels
end
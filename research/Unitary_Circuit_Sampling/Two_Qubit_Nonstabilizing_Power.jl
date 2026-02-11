"""
Two_Qubit_Nonstabilizing_Power.jl

Estimate the "non-stabilising power" of random 2-qubit unitaries by applying them
to a set of stabiliser states and measuring the resulting magic.

File output
- Disabled by default.
- Enable with: SAVE_OUTPUT=1 julia Two_Qubit_Nonstabilizing_Power.jl
- Outputs go to research/output/Unitary_Circuit_Sampling/Two_Qubit_Nonstabilizing_Power/

Notes
- This script still relies on legacy Modules/*.jl includes. We keep the logic intact here.
"""

using Random
using ProgressBars
using JLD2
using Base.Threads
using MAT

include(joinpath(@__DIR__, "..", "research_utils.jl"))

current_dir = @__DIR__

# Legacy modules.
filepath = joinpath(current_dir, "..", "Modules", "Magic.jl")
include(filepath)
filepath = joinpath(current_dir, "..", "Modules", "Random_Unitaries.jl")
include(filepath)
filepath = joinpath(current_dir, "..", "Modules", "Entanglement.jl")
include(filepath)

using .Measure_Entanglement
using .Measure_Magic
using .Random_Unitary_Generation

function main()
    save_output = get_bool_env("SAVE_OUTPUT", false)

    seed = 1
    Random.seed!(seed)

    no_samples = 2^21
    no_qubits = 2

    # Generate all 2-qubit Clifford gates (legacy helper from this script).
    cliffords, labels = Generate_All_2_Qubit_Clifford_Gates()

    zero_state = [1.0 + 0im, 0.0 + 0im, 0.0 + 0im, 0.0 + 0im]
    clifford_states = [c * zero_state for c in cliffords]

    function unique_states(states)
        unique_vecs = []
        for v in states
            is_unique = true
            for u in unique_vecs
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

    psi_list = unique_states([round.(v; digits=10) for v in clifford_states])
    no_stab = length(psi_list)

    psi_mat = hcat(psi_list...)'  # each state is a row

    strings = Measure_Magic.GenerateAllPauliStrings(2)
    pauli_ops = Measure_Magic.PauliOperatorList(strings, 2)

    magic = zeros(no_samples, no_stab)

    for i in ProgressBar(1:no_samples)
        u = Random_Unitary_Generation.Generate_Regular_Unitary_Circuit(no_qubits)
        for stabstate in 1:no_stab
            psi_0 = psi_mat[stabstate, :]
            state = u * psi_0
            magic[i, stabstate] = Measure_Magic.MeasureMagic_Pure(state, pauli_ops, 2)[1]
        end
    end

    if save_output
        out_mat = output_path(
            "NonStabilisingPower_of_Unitary_N_2_Samples_$(no_samples)_Seed_$(seed).mat";
            script_dir=@__DIR__,
            script_file=@__FILE__,
        )
        ensure_parent_dir(out_mat)
        MAT.matwrite(out_mat, Dict("Magic" => magic))

        out_jld2 = output_path(
            "UnitaryMagicSampling_AllStabs_N_2_Samples_$(no_samples).jld2";
            script_dir=@__DIR__,
            script_file=@__FILE__,
        )
        ensure_parent_dir(out_jld2)
        JLD2.@save out_jld2 magic

        @info "Saved outputs to $out_mat and $out_jld2"
    else
        @info "SAVE_OUTPUT is disabled; not writing MAT/JLD2 files."
    end

    return nothing
end

# Original helper copied verbatim (kept outside main for clarity)
function Generate_All_2_Qubit_Clifford_Gates()
    Id = [1 0; 0 1]
    PX = [0 1; 1 0]
    PY = [0 -1im; 1im 0]
    PZ = [1 0; 0 -1]
    H = 1 / sqrt(2) * [1 1; 1 -1]
    S = [1 0; 0 1im]

    clifford_gates = []

    had = [Id, H]
    had_labels = ["I", "H"]

    phase = [Id, H * S, S * H]
    phase_labels = ["I", "HS", "SH"]

    pauli = [Id, PX, PY, PZ]
    pauli_labels = ["I", "X", "Y", "Z"]

    labels = String[]

    for hj1 in 1:length(had), hj2 in 1:length(had),
        vj1 in 1:length(phase), vj2 in 1:length(phase),
        pj1 in 1:length(pauli), pj2 in 1:length(pauli)

        push!(clifford_gates, kron(had[hj1], had[hj2]) * kron(phase[vj1], phase[vj2]) * kron(pauli[pj1], pauli[pj2]))
        push!(labels, "($(had_labels[hj1])⊗$(had_labels[hj2]))*" *
                      "($(phase_labels[vj1])⊗$(phase_labels[vj2]))*" *
                      "($(pauli_labels[pj1])⊗$(pauli_labels[pj2]))")
    end

    cnot = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]

    for hj1 in 1:length(had), hj2 in 1:length(had),
        vj1 in 1:length(phase), vj2 in 1:length(phase),
        vj3 in 1:length(phase), vj4 in 1:length(phase),
        pj1 in 1:length(pauli), pj2 in 1:length(pauli)

        push!(clifford_gates, kron(had[hj1], had[hj2]) * kron(phase[vj1], phase[vj2]) * cnot * kron(phase[vj3], phase[vj4]) * kron(pauli[pj1], pauli[pj2]))
        push!(labels, "($(had_labels[hj1])⊗$(had_labels[hj2]))*" *
                      "($(phase_labels[vj1])⊗$(phase_labels[vj2]))*" *
                      "CNOT*" *
                      "($(phase_labels[vj3])⊗$(phase_labels[vj4]))*" *
                      "($(pauli_labels[pj1])⊗$(pauli_labels[pj2]))")
    end

    swap = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]

    for hj1 in 1:length(had), hj2 in 1:length(had),
        vj1 in 1:length(phase), vj2 in 1:length(phase),
        vj3 in 1:length(phase), vj4 in 1:length(phase),
        pj1 in 1:length(pauli), pj2 in 1:length(pauli)

        push!(clifford_gates, kron(had[hj1], had[hj2]) * kron(phase[vj1], phase[vj2]) * swap * cnot * kron(phase[vj3], phase[vj4]) * kron(pauli[pj1], pauli[pj2]))
        push!(labels, "($(had_labels[hj1])⊗$(had_labels[hj2]))*" *
                      "($(phase_labels[vj1])⊗$(phase_labels[vj2]))*" *
                      "SWAP*CNOT*" *
                      "($(phase_labels[vj3])⊗$(phase_labels[vj4]))*" *
                      "($(pauli_labels[pj1])⊗$(pauli_labels[pj2]))")
    end

    for hj1 in 1:length(had), hj2 in 1:length(had),
        vj1 in 1:length(phase), vj2 in 1:length(phase),
        pj1 in 1:length(pauli), pj2 in 1:length(pauli)

        push!(clifford_gates, kron(had[hj1], had[hj2]) * kron(phase[vj1], phase[vj2]) * swap * kron(pauli[pj1], pauli[pj2]))
        push!(labels, "($(had_labels[hj1])⊗$(had_labels[hj2]))*" *
                      "($(phase_labels[vj1])⊗$(phase_labels[vj2]))*" *
                      "SWAP*" *
                      "($(pauli_labels[pj1])⊗$(pauli_labels[pj2]))")
    end

    return clifford_gates, labels
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

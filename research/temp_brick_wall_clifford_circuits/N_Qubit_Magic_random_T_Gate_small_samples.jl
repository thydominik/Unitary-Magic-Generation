"""
N_Qubit_Magic_random_T_Gate_small_samples.jl

Same basic experiment as N_Qubit_Magic_random_T_Gate.jl, but with smaller sample counts
for larger N.

File output
- Disabled by default.
- Enable with: SAVE_OUTPUT=1 julia N_Qubit_Magic_random_T_Gate_small_samples.jl

Outputs go to research/output/TEMP_BrickWall_Clifford/N_Qubit_Magic_random_T_Gate_small_samples/
"""

using LinearAlgebra
using SparseArrays
using Random
using ProgressBars
using JLD2
using Base.Threads
using MAT

include(joinpath(@__DIR__, "..", "research_utils.jl"))

# The helper functions are identical to the main script; for now we keep them in-place.
# If you want, we can de-duplicate later by moving them into a shared TEMP helper.

include(joinpath(@__DIR__, "N_Qubit_Magic_random_T_Gate.jl"))

function main()
    save_output = get_bool_env("SAVE_OUTPUT", false)

    # The included file defines the helper functions; we reuse them here.
    current_dir = @__DIR__

    include(joinpath(current_dir, "..", "Modules", "Magic.jl"))
    include(joinpath(current_dir, "..", "Modules", "Random_Unitaries.jl"))

    using .Random_Unitary_Generation
    using .Measure_Magic

    for NoQ in 2:8
        for Depth in (5 * NoQ,)
            no_qubits = NoQ
            d = Depth

            clf_gates, _ = Generate_All_2_Qubit_Clifford_Gates()
            strings = Measure_Magic.GenerateAllPauliStrings(no_qubits)
            pauli_ops = Measure_Magic.PauliOperatorList(strings, no_qubits)

            Random.seed!(1)
            no_samples = 2^12
            psi_0 = zeros(2^no_qubits)
            psi_0[1] = 1

            for nT in 1:d
                magic_vals = Vector{Float64}()

                for _ in ProgressBar(1:no_samples)
                    possible = generate_T_gate_coordinates(no_qubits, d)
                    placement = zeros(Int, nT, 2)

                    for i in 1:nT
                        idx = rand(1:size(possible, 1))
                        placement[i, :] = possible[idx, :]
                        possible = delete_row(possible, idx)
                    end

                    gate = Build_CLF_BrickWall_Circuit_with_T(no_qubits, d, placement, clf_gates)
                    psi = gate * psi_0
                    push!(magic_vals, Measure_Magic.MeasureMagic_Pure(psi, pauli_ops, 2))
                end

                if save_output
                    fname = "CLF_BW_N_$(no_qubits)_D_$(d)_NT_$(nT)_Samples_$(no_samples).mat"
                    out_mat = output_path(fname; script_dir=@__DIR__, script_file=@__FILE__)
                    ensure_parent_dir(out_mat)
                    MAT.matwrite(out_mat, Dict("Magic" => magic_vals))
                end
            end
        end
    end

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

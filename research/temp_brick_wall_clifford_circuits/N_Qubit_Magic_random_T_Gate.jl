"""
N_Qubit_Magic_random_T_Gate.jl

Brickwall-style circuit experiment with random 2-qubit Clifford gates and inserted T gates.

This folder is marked TEMP, so the code is kept close to the original.
The main changes are:
- Safe-by-default output: no files are written unless SAVE_OUTPUT=1.
- Outputs go under research/output/TEMP_BrickWall_Clifford/N_Qubit_Magic_random_T_Gate/

Notes
- This script currently uses legacy Modules/*.jl includes.
- It can be very slow for large N, depth, and sample counts.
"""

using LinearAlgebra
using SparseArrays
using Random
using ProgressBars
using JLD2
using Base.Threads
using MAT

include(joinpath(@__DIR__, "..", "research_utils.jl"))

# --- Original helper functions (kept verbatim) ---------------------------------
# (Generate_All_2_Qubit_Clifford_Gates, generate_T_gate_coordinates, etc.)

function Generate_All_2_Qubit_Clifford_Gates()
    Id = sparse([1 0; 0 1])
    PX = sparse([0 1; 1 0])
    PY = sparse([0 -im; im 0])
    PZ = sparse([1 0; 0 -1])
    H = 1 / sqrt(2) * [1 1; 1 -1]
    S = sparse([1 0; 0 im])

    Clifford_Gates_2_Site = []

    Hadamard_group = [Id, H]
    Hadamard_group_labels = ["I", "H"]

    Phaseshift_group = [Id, H * S, S * H]
    Phaseshift_group_labels = ["I", "HS", "SH"]

    Pauli_group = [Id, PX, PY, PZ]
    Pauli_group_labels = ["I", "X", "Y", "Z"]

    Dimension_Hadamard_group = length(Hadamard_group)
    Dimension_Phaseshift_group = length(Phaseshift_group)
    Dimension_Pauli_group = length(Pauli_group)

    Clifford_Lables = String[]
    for hj1 in 1:Dimension_Hadamard_group, hj2 in 1:Dimension_Hadamard_group,
        vj1 in 1:Dimension_Phaseshift_group, vj2 in 1:Dimension_Phaseshift_group,
        pj1 in 1:Dimension_Pauli_group, pj2 in 1:Dimension_Pauli_group

        push!(Clifford_Gates_2_Site, sparse(kron(Hadamard_group[hj1], Hadamard_group[hj2]) *
                                           kron(Phaseshift_group[vj1], Phaseshift_group[vj2]) *
                                           kron(Pauli_group[pj1], Pauli_group[pj2])))
        push!(Clifford_Lables, "($(Hadamard_group_labels[hj1])⊗$(Hadamard_group_labels[hj2]))*" *
                               "($(Phaseshift_group_labels[vj1])⊗$(Phaseshift_group_labels[vj2]))*" *
                               "($(Pauli_group_labels[pj1])⊗$(Pauli_group_labels[pj2]))")
    end

    CNOT = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]

    for hj1 in 1:Dimension_Hadamard_group, hj2 in 1:Dimension_Hadamard_group,
        vj1 in 1:Dimension_Phaseshift_group, vj2 in 1:Dimension_Phaseshift_group,
        vj3 in 1:Dimension_Phaseshift_group, vj4 in 1:Dimension_Phaseshift_group,
        pj1 in 1:Dimension_Pauli_group, pj2 in 1:Dimension_Pauli_group

        push!(Clifford_Gates_2_Site, sparse(kron(Hadamard_group[hj1], Hadamard_group[hj2]) *
                                           kron(Phaseshift_group[vj1], Phaseshift_group[vj2]) *
                                           CNOT *
                                           kron(Phaseshift_group[vj3], Phaseshift_group[vj4]) *
                                           kron(Pauli_group[pj1], Pauli_group[pj2])))
        push!(Clifford_Lables, "($(Hadamard_group_labels[hj1])⊗$(Hadamard_group_labels[hj2]))*" *
                               "($(Phaseshift_group_labels[vj1])⊗$(Phaseshift_group_labels[vj2]))*" *
                               "CNOT*" *
                               "($(Phaseshift_group_labels[vj3])⊗$(Phaseshift_group_labels[vj4]))*" *
                               "($(Pauli_group_labels[pj1])⊗$(Pauli_group_labels[pj2]))")
    end

    SWAP = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]

    for hj1 in 1:Dimension_Hadamard_group, hj2 in 1:Dimension_Hadamard_group,
        vj1 in 1:Dimension_Phaseshift_group, vj2 in 1:Dimension_Phaseshift_group,
        vj3 in 1:Dimension_Phaseshift_group, vj4 in 1:Dimension_Phaseshift_group,
        pj1 in 1:Dimension_Pauli_group, pj2 in 1:Dimension_Pauli_group

        push!(Clifford_Gates_2_Site, sparse(kron(Hadamard_group[hj1], Hadamard_group[hj2]) *
                                           kron(Phaseshift_group[vj1], Phaseshift_group[vj2]) *
                                           SWAP * CNOT *
                                           kron(Phaseshift_group[vj3], Phaseshift_group[vj4]) *
                                           kron(Pauli_group[pj1], Pauli_group[pj2])))
        push!(Clifford_Lables, "($(Hadamard_group_labels[hj1])⊗$(Hadamard_group_labels[hj2]))*" *
                               "($(Phaseshift_group_labels[vj1])⊗$(Phaseshift_group_labels[vj2]))*" *
                               "SWAP*CNOT*" *
                               "($(Phaseshift_group_labels[vj3])⊗$(Phaseshift_group_labels[vj4]))*" *
                               "($(Pauli_group_labels[pj1])⊗$(Pauli_group_labels[pj2]))")
    end

    for hj1 in 1:Dimension_Hadamard_group, hj2 in 1:Dimension_Hadamard_group,
        vj1 in 1:Dimension_Phaseshift_group, vj2 in 1:Dimension_Phaseshift_group,
        pj1 in 1:Dimension_Pauli_group, pj2 in 1:Dimension_Pauli_group

        push!(Clifford_Gates_2_Site, sparse(kron(Hadamard_group[hj1], Hadamard_group[hj2]) *
                                           kron(Phaseshift_group[vj1], Phaseshift_group[vj2]) *
                                           SWAP *
                                           kron(Pauli_group[pj1], Pauli_group[pj2])))
        push!(Clifford_Lables, "($(Hadamard_group_labels[hj1])⊗$(Hadamard_group_labels[hj2]))*" *
                               "($(Phaseshift_group_labels[vj1])⊗$(Phaseshift_group_labels[vj2]))*" *
                               "SWAP*" *
                               "($(Pauli_group_labels[pj1])⊗$(Pauli_group_labels[pj2]))")
    end

    return Clifford_Gates_2_Site, Clifford_Lables
end

function generate_T_gate_coordinates(N::Int, D::Int)
    coordinates = zeros(Int, (N - 1) * D, 2)
    k = 1

    for d in 1:D
        if isodd(d)
            for n in 1:N
                if n == N && iseven(n)
                    coordinates[k, :] = [n, d]
                    k += 1
                elseif n < N
                    coordinates[k, :] = [n, d]
                    k += 1
                end
            end
        else
            for n in 2:N
                if n == N && isodd(n)
                    coordinates[k, :] = [n, d]
                    k += 1
                elseif n < N
                    coordinates[k, :] = [n, d]
                    k += 1
                end
            end
        end
    end

    return coordinates
end

function Build_CLF_BrickWall_Circuit_with_T(N::Int, D::Int, T_coords, CLFGates)
    TGate = sparse([[1 0]; [0 exp(im * pi / 4)]])
    Gate = sparse(zeros(Float64, 2^N, 2^N))

    for depthIndex in 1:D
        if mod(depthIndex, 2) == 1
            Layer = sparse(I(4))
            for qubitIndex in 1:2:N
                if qubitIndex < N
                    RandomGate = CLFGates[rand(1:11520)]
                    if qubitIndex == 1
                        Layer = Layer * RandomGate
                    else
                        Layer = kron(Layer, RandomGate)
                    end
                else
                    Layer = kron(Layer, I(2))
                end
            end
        else
            Layer = I(2)
            for qubitIndex in 2:2:N
                if qubitIndex < N
                    RandomGate = CLFGates[rand(1:11520)]
                    Layer = kron(Layer, RandomGate)
                else
                    Layer = kron(Layer, I(2))
                end
            end
        end

        Gate = depthIndex == 1 ? Layer : Layer * Gate

        if sum(T_coords[:, 2] .== depthIndex) > 0
            for i in 1:size(T_coords, 1)
                if T_coords[i, 2] == depthIndex
                    extra = kron(kron(I(2^(T_coords[i, 1] - 1)), TGate), I(2^(N - T_coords[i, 1])))
                    Gate = extra * Gate
                end
            end
        end
    end

    return Gate
end

function delete_row(matrix, row_index)
    return vcat(matrix[1:row_index-1, :], matrix[row_index+1:end, :])
end

# --- Main experiment -----------------------------------------------------------

function main()
    save_output = get_bool_env("SAVE_OUTPUT", false)

    current_dir = @__DIR__

    # Legacy includes.
    include(joinpath(current_dir, "..", "Modules", "Magic.jl"))
    include(joinpath(current_dir, "..", "Modules", "Random_Unitaries.jl"))

    using .Random_Unitary_Generation
    using .Measure_Magic

    for NoQ = 2:8
        for Depth in (5 * NoQ, 10 * NoQ)
            no_qubits = NoQ
            d = Depth

            clf_gates, _ = Generate_All_2_Qubit_Clifford_Gates()
            strings = Measure_Magic.GenerateAllPauliStrings(no_qubits)
            pauli_ops = Measure_Magic.PauliOperatorList(strings, no_qubits)

            Random.seed!(1)
            no_samples = 2^16
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

    @info "Done. (Set SAVE_OUTPUT=1 if you want MAT files written.)"
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

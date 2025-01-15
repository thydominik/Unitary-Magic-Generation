using LinearAlgebra
using SparseArrays
using Random
using ProgressBars
using JLD2
using Base.Threads
using MAT

function Generate_All_2_Qubit_Clifford_Gates()
    # Identity
    Id = sparse([1 0; 0 1])
    # Pauli Gates
    PX = sparse([0 1; 1 0])
    PY = sparse([0 -im; im 0])
    PZ = sparse([1 0; 0 -1])
    # Hadamard
    H = 1/sqrt(2) * [1 1; 1 -1]
    # Phase gate
    S = sparse([1 0; 0 im])
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

    # these will be accessed multiple times
    Dimension_Hadamard_group    = length(Hadamard_group)
    Dimension_Phaseshift_group  = length(Phaseshift_group)
    Dimension_Pauli_group       = length(Pauli_group)

    # Gates created from the 1 qubit clifford gates (2, 2, 3, 3, 4, 4)
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

    # CNOT class contains 5184 elements = 24^2*3^2 or (2, 2, 3, 3, 3, 3, 4, 4)
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

    # class 3 elements with SWAP and CNOT class contains 5184 elements
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

    # class 4 elements contains 576 elements (2, 2, 3, 3, 4, 4).
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

    Clifford_Gates_2_Qubits = Clifford_Gates_2_Site
    Labels = Clifford_Lables
    return Clifford_Gates_2_Qubits, Labels
end

function generate_T_gate_coordinates(N::Int, D::Int)
    #coordinates = Tuple{Int, Int}[]
    coordinates = zeros(Int, ( N - 1)*D, 2)
    k = 1
    for d in 1:D
        if isodd(d)
            for n in 1:N
                if n == N && iseven(n)
                    coordinates[k, :] = [n, d];
                    k += 1
                elseif n<N
                    coordinates[k, :] = [n, d];
                    k += 1
                end
            end
        else
            for n in 2:N
                if n == N && isodd(n)
                    coordinates[k, :] = [n, d];
                    k += 1
                elseif n<N
                    coordinates[k, :] = [n, d];
                    k += 1
                end
            end
        end
    end
    
    return coordinates
end

function Build_CLF_BrickWall_Circuit_with_T(N::Int, D::Int, T_coords, CLFGates)

    TGate = sparse([[1 0]; [0 exp(im * pi/4)]]);
    Gate = sparse(zeros(Float64, 2^N, 2^N))
    for depthIndex in 1:D
        if mod(depthIndex, 2) == 1
            Layer = sparse(I(4))

            for qubitIndex in 1:2:N
                if qubitIndex < N
                    RandomGate = CLFGates[rand(1:11520)];
                    if qubitIndex == 1
                        Layer = Layer * RandomGate;
                    else
                        Layer = kron(Layer, RandomGate);
                    end
                else
                    RandomGate = I(2);
                    Layer = kron(Layer, RandomGate);
                end
            end
        else
            Layer = I(2);
            for qubitIndex in 2:2:N
                if qubitIndex < N
                    RandomGate = CLFGates[rand(1:11520)];
                    Layer = kron(Layer, RandomGate);
                elseif qubitIndex == N
                    RandomGate = I(2);
                    Layer = kron(Layer, RandomGate);
                end
            end
        end

        if depthIndex == 1
            Gate = Layer;
        else
            Gate = Layer * Gate;
        end

        if sum(T_coords[:, 2] .== depthIndex) > 0
            for i in 1:size(T_coords, 1)
                if T_coords[i, 2] == depthIndex
                    ExtraLayer = kron(kron(I(2^(T_coords[i, 1] - 1)), TGate), I(2^(N - T_coords[i, 1])));
                    Gate = ExtraLayer * Gate;
                end
            end
        end
    end
    return Gate
end


function delete_row(matrix, row_index)
    return vcat(matrix[1:row_index-1, :], matrix[row_index+1:end, :])
end

TGate = sparse([[1, 0], [0, exp(im * pi/4)]])

current_dir = @__DIR__

# Unitary matrices
filepath = joinpath(current_dir, "..", "Modules", "Magic.jl")
include(filepath)
# Magic
filepath = joinpath(current_dir, "..", "Modules", "Random_Unitaries.jl")
include(filepath)
#include("C:\\Dominik\\PhD\\Projects\\Unitary-Magic-Generation\\Modules\\Magic.jl")
using .Random_Unitary_Generation
using .Measure_Magic

for NoQ = 2:9
    for Depth = [5*NoQ, 10*NoQ]
        No_Qubits = NoQ;
        D = Depth;

        CLFGates, CLFLabels = Generate_All_2_Qubit_Clifford_Gates();
        TGate = sparse([[1 0]; [0 exp(im * pi/4)]]);
        Strings = Measure_Magic.GenerateAllPauliStrings(No_Qubits);
        PauliOperators = Measure_Magic.PauliOperatorList(Strings, No_Qubits);
        Random.seed!(1);
        No_Samples = 2^16
        #Psi_0 = 1/sqrt(2^No_Qubits) * ones(2^No_Qubits);
        Psi_0 = zeros(2^No_Qubits); Psi_0[1] = 1;

        for nT in 1:D
            Magic = Vector{Float64}()
            for iterations in ProgressBar(1:No_Samples)
                PossibleTGateCoordinates = generate_T_gate_coordinates(No_Qubits, D);
                TGatePlacement = zeros(Int, nT, 2)
                for i in 1:nT
                    TGateIndex = rand(1:size(PossibleTGateCoordinates, 1));
                    TGatePlacement[i, :] = PossibleTGateCoordinates[TGateIndex, :];
                    PossibleTGateCoordinates = delete_row(PossibleTGateCoordinates, TGateIndex);
                end

                Gate = Build_CLF_BrickWall_Circuit_with_T(No_Qubits, D, TGatePlacement, CLFGates)
                Psi = Gate * Psi_0
                push!(Magic, Measure_Magic.MeasureMagic_Pure(Psi, PauliOperators, 2))
            end
            fname = "CLF_BW_N_$(No_Qubits)_D_$(D)_NT_$(nT)_Samples_$(No_Samples).mat"
            matwrite(fname, Dict("Magic" => Magic))

        end
    end
end
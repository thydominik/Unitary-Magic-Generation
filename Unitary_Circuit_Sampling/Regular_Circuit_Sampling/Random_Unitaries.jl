

module Random_Unitary_Generation

using LinearAlgebra
using RandomMatrices

export CUE_Matrix, Generate_BW_Unitary_Circuit, Generate_Regular_Unitary_Circuit

function CUE_Matrix(Matrix_Size::Int, beta = 2)
    U = rand(Haar(beta), Matrix_Size)
    return U
end

function Generate_Regular_Unitary_Circuit(No_Qubits::Int)
    Operator = CUE_Matrix(2^No_Qubits, 2)
end

function Generate_BW_Unitary_Circuit(No_Qubits::Int, Depth::Int = 5 * No_Qubits)
    Circuit = [];
    for DepthIndex in 1:Depth
        if DepthIndex % 2 == 1
            Layer = I(4);
            for QubitIndex in 1:2:No_Qubits
                if QubitIndex < No_Qubits
                    RandomUnitary = CUE_Matrix(4);
                    if QubitIndex == 1
                        Layer = Layer * RandomUnitary;
                    else
                        Layer = kron(Layer, RandomUnitary);
                    end # IF qubitIndex == 1
                else # IF qubitIndex < No_Qubits
                    RandomUnitary = I(2);
                    Layer = kron(Layer, RandomUnitary);
                end # IF qubitIndex < No_Qubits
            end # FOR qubitIndex
        else    # IF depthIndex % 2
            Layer = I(2);
            for qubitIndex in 2:2:No_Qubits
                if qubitIndex < No_Qubits
                    RandomUnitary = CUE_Matrix(4);
                    Layer = kron(Layer, RandomUnitary);
                elseif qubitIndex == No_Qubits
                    RandomUnitary = I(2);
                    Layer = kron(Layer, RandomUnitary);
                end # IF qubitIndex < No_Qubits 
            end # FOR qubitIndex
        end # IF depthIndex % 2
        if DepthIndex == 1
            Circuit = Layer
        else
            Circuit = Layer * Circuit
        end
    end # FOR DepthIndex

    return Circuit
end
end
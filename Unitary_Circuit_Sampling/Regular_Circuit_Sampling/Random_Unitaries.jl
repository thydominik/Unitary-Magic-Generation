using RandomMatrices

module Random_Unitary_Generation

function CUE_Matrix(Matrix_Size, beta)
    U = rand(Haar(beta), Matrix_Size)
    return U
end

function Generate_Regular_Unitary_Circuit(No_Qubits)
    Operator = CUE_Matrix(Matrix_Size = 2^N, beta = 2)
end

end
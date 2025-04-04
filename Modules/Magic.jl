module Measure_Magic

using LinearAlgebra
using SparseArrays

export GenerateAllPauliStrings, PauliMatrix, MeasureMagic_Pure, MeasureMagic_Mixed, PauliOperatorList

function PauliOperatorList(Strings, No_Qubits::Int)
    σ = Vector{SparseMatrixCSC{ComplexF64}}(undef, 4^No_Qubits)
    for PauliIndex in 1:4^No_Qubits
        CurrentString   = Strings[PauliIndex]
        Operator        = PauliMatrix(CurrentString[1])
        for QubitIndex in 2:No_Qubits
            LocalPauliMatrix    = PauliMatrix(CurrentString[QubitIndex])
            Operator            = kron(Operator, LocalPauliMatrix)
        end # FOR QubitIndex
        σ[PauliIndex] = sparse(Operator)
    end # FOR PauliIndex
    return σ
end

function MeasureMagic_Mixed(ϱ, σ, α = 2)
    # Calculating the stabiliser Rényi Entropy
    # https://doi.org/10.1103/PhysRevLett.128.050402

    No_Qubits = Int(log2(size(ϱ, 1)))

    HilbertSpaceDimension = 2^No_Qubits
    
    Ξ = Vector{Float64}(undef, 4^No_Qubits)
    for PauliIndex in 1:4^No_Qubits
        Operator = σ[PauliIndex]
        Ξ[PauliIndex] = real((1/HilbertSpaceDimension) * sum(diag(Operator * ϱ))^2)
    end # FOR PauliIndex
    Magic = (1 - α)^(-1) * log2(sum(Ξ.^α)) - log2(HilbertSpaceDimension) + log2(sum(eigvals(ϱ).^2))
    return Magic
end

function MeasureMagic_Pure(State, PauliOperators, α::Union{Integer, Float64, Vector{<:Number}} = 2)
    No_Qubits = Int(log2(length(State)))
    HilbertSpaceDimension = 2^No_Qubits
    
    Ξ = Vector{Float64}(undef, 4^No_Qubits)
    for PauliIndex in 1:4^No_Qubits
        Operator = PauliOperators[PauliIndex]
        Ξ[PauliIndex] = real((1/HilbertSpaceDimension) * (conj(transpose(State[:])) * Operator * State[:])^2)
    end
    
    if α isa Number
        if α == 1
            Magic = 1 - HilbertSpaceDimension * sum(Ξ.^2)
        else
            Magic = (1 - float(α))^(-1) * log2(sum(Ξ.^float(α))) - log2(HilbertSpaceDimension)
        end
        
    else
        Magic = [(1 - float(α_i))^(-1) * log2(sum(Ξ.^float(α_i))) - log2(HilbertSpaceDimension) for α_i in α]
    end
    
    return Magic, Ξ
end

function GenerateAllPauliStrings(No_Qubits::Int)
    Pauli_Matrices = ["I", "X", "Y", "Z"]
    Pauli_Strings = [join(string) for string in Iterators.product(fill(Pauli_Matrices, No_Qubits)...)]
    return Pauli_Strings
end

function PauliMatrix(which_Pauli)
    if which_Pauli == 'I' || which_Pauli == 0
        Pauli_Matrix = sparse([1.0 0.0;0.0 1.0])
    elseif which_Pauli == 'X' || which_Pauli == 1
        Pauli_Matrix = sparse([0.0 1.0; 1.0 0.0])
    elseif which_Pauli == 'Y' || which_Pauli == 2
        Pauli_Matrix = sparse([0.0 -im; im 0.0])
    elseif which_Pauli == 'Z' || which_Pauli == 3
        Pauli_Matrix = sparse([1.0 0.0; 0.0 -1.0])
    end

    return Pauli_Matrix
end

end
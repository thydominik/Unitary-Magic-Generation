module Measure_Magic

using LinearAlgebra

export GenerateAllPauliStrings, PauliMatrix, MeasureMagic

function MeasureMagic(State, Strings, α = 2)
    # Calculating the stabiliser Rényi Entropy
    # https://doi.org/10.1103/PhysRevLett.128.050402

    No_Qubits = Int(log2(length(State)))

    HilbertSpaceDimension = 2^No_Qubits
    Ξ = Vector()
    for PauliIndex in 1:4^No_Qubits
        CurrentString   = Strings[PauliIndex]
        Operator        = PauliMatrix(CurrentString[1])
        for QubitIndex in 2:No_Qubits
            LocalPauliMatrix    = PauliMatrix(CurrentString[QubitIndex])
            Operator            = kron(Operator, LocalPauliMatrix)
        end # FOR QubitIndex
        push!(Ξ, (1/HilbertSpaceDimension) * (conj(transpose(State[:])) * Operator * State[:])^2)
    end # FOR PauliIndex
    Magic = real((1 - α)^(-1) * log(sum(Ξ.^α)) - log(HilbertSpaceDimension))
    return Magic
end

function GenerateAllPauliStrings(No_Qubits::Int)
    Pauli_Matrices = ["I", "X", "Y", "Z"]
    Pauli_Strings = [join(string) for string in Iterators.product(fill(Pauli_Matrices, No_Qubits)...)]
    return Pauli_Strings
end

function PauliMatrix(which_Pauli)
    if which_Pauli == 'I' || which_Pauli == 0
        Pauli_Matrix = [1.0 0.0;0.0 1.0]
    elseif which_Pauli == 'X' || which_Pauli == 1
        Pauli_Matrix = [0.0 1.0; 1.0 0.0]
    elseif which_Pauli == 'Y' || which_Pauli == 2
        Pauli_Matrix = [0.0 -im; im 0.0]
    elseif which_Pauli == 'Z' || which_Pauli == 3
        Pauli_Matrix = [1.0 0.0; 0.0 -1.0]
    end

    return Pauli_Matrix
end

end
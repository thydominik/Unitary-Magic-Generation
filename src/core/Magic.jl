"""
    Magic

Core routines for magic quantification and Pauli operator handling.

This module provides utilities to construct Pauli operator strings and matrices and to
compute magic measures for pure and mixed quantum states based on stabiliser Rényi entropy.

The implementation closely follows:
- https://doi.org/10.1103/PhysRevLett.128.050402

# Exported API

- `pauli_operator_list` — build tensor-product Pauli operators from string labels
- `pauli_matrix`       — single-qubit Pauli matrix (sparse)
- `measure_magic_pure` — magic of pure states
- `measure_magic_mixed` — magic of mixed states
- `generate_all_pauli_strings` — generate all Pauli strings for N qubits
"""
module Magic

using LinearAlgebra
using SparseArrays

export 
    pauli_operator_list,
    pauli_matrix,
    measure_magic_pure,
    measure_magic_mixed,
    generate_all_pauli_strings

"""
    pauli_operator_list(strings::Vector{<:AbstractString}, n_qubits::Integer) -> Vector{SparseMatrixCSC{ComplexF64,Int}}

Construct the list of tensor-product Pauli operators corresponding to the given
Pauli strings.

Each string is interpreted as a sequence of characters from the set `{"I","X","Y","Z"}`
for each qubit.

# Arguments
- `strings`: Vector of Pauli strings of length `n_qubits` each
- `n_qubits`: Number of qubits

# Returns
- Vector of sparse matrices of size `2^n_qubits × 2^n_qubits`
"""
function pauli_operator_list(strings::Vector{<:AbstractString}, n_qubits::Integer)
    n_qubits ≤ 0 && throw(ArgumentError("n_qubits must be positive, got $n_qubits"))

    σ = Vector{SparseMatrixCSC{ComplexF64,Int}}(undef, length(strings))
    for (i, current_string) in enumerate(strings)
        length(current_string) == n_qubits || throw(ArgumentError(
            "Pauli string length $(length(current_string)) does not match n_qubits=$n_qubits",
        ))

        op = pauli_matrix(current_string[1])
        for q in 2:n_qubits
            local_pauli = pauli_matrix(current_string[q])
            op = kron(op, local_pauli)
        end
        σ[i] = sparse(op)
    end
    return σ
end

"""
    measure_magic_mixed(ρ::AbstractMatrix{<:Complex}, σ::Vector{<:SparseMatrixCSC}; α::Real=2) -> Float64

Compute the magic of a mixed state ρ using stabiliser Rényi entropy.

# Arguments
- `ρ`: Density matrix of size `2^n × 2^n`
- `σ`: Vector of Pauli operators (see `pauli_operator_list`)
- `α`: Rényi parameter (default `2`)

# Returns
- Magic value as `Float64`
"""
function measure_magic_mixed(ρ::AbstractMatrix{<:Complex}, σ::Vector{<:SparseMatrixCSC}; α::Real=2)
    n_qubits = Int(log2(size(ρ, 1)))
    hilbert_dim = 2^n_qubits

    Ξ = Vector{Float64}(undef, length(σ))
    for (i, op) in enumerate(σ)
        Ξ[i] = real((1 / hilbert_dim) * sum(diag(op * ρ))^2)
    end

    magic = (1 - α)^(-1) * log2(sum(Ξ .^ α)) - log2(hilbert_dim) + log2(sum(eigvals(ρ) .^ 2))
    return magic
end

"""
    measure_magic_pure(state::AbstractVector{<:Complex}, pauli_ops::Vector{<:SparseMatrixCSC}, α=2) -> (magic, Ξ)

Compute the magic of a pure state using stabiliser Rényi entropy.

# Arguments
- `state`: State vector of length `2^n`
- `pauli_ops`: Vector of Pauli operators
- `α`: Can be a single number or a vector of parameters

# Returns
- `magic`: Scalar or vector of magic values depending on `α`
- `Ξ`: Vector of stabiliser probabilities
"""
function measure_magic_pure(
    state::AbstractVector{<:Complex},
    pauli_ops::Vector{<:SparseMatrixCSC},
    α::Union{Integer, Float64, AbstractVector{<:Number}}=2,
)
    n_qubits = Int(log2(length(state)))
    hilbert_dim = 2^n_qubits

    Ξ = Vector{Float64}(undef, length(pauli_ops))
    ψ = state[:]
    ψ† = conj(transpose(ψ))

    for (i, op) in enumerate(pauli_ops)
        Ξ[i] = real((1 / hilbert_dim) * (ψ† * op * ψ)^2)
    end

    if α isa Number
        if α == 1
            magic = 1 - hilbert_dim * sum(Ξ .^ 2)
        else
            a = float(α)
            magic = (1 - a)^(-1) * log2(sum(Ξ .^ a)) - log2(hilbert_dim)
        end
    else
        magic = [(1 - float(a))^(-1) * log2(sum(Ξ .^ float(a))) - log2(hilbert_dim) for a in α]
    end

    return magic, Ξ
end

"""
    generate_all_pauli_strings(n_qubits::Integer) -> Vector{String}

Generate all Pauli strings for `n_qubits` qubits using the alphabet `{"I","X","Y","Z"}`.

# Arguments
- `n_qubits`: Number of qubits

# Returns
- Vector of strings of length `n_qubits`, each a tensor-product Pauli label
"""
function generate_all_pauli_strings(n_qubits::Integer)
    n_qubits ≤ 0 && throw(ArgumentError("n_qubits must be positive, got $n_qubits"))
    pauli_letters = ["I", "X", "Y", "Z"]
    return [join(s) for s in Iterators.product(fill(pauli_letters, n_qubits)...)]
end

"""
    pauli_matrix(which::Union{Char,Integer}) -> SparseMatrixCSC{ComplexF64,Int}

Return the single-qubit Pauli matrix corresponding to the given label.

Accepted labels:
- `'I'` or `0`
- `'X'` or `1`
- `'Y'` or `2`
- `'Z'` or `3`
"""
function pauli_matrix(which::Union{Char,Integer})::SparseMatrixCSC{ComplexF64,Int}
    if which == 'I' || which == 0
        return sparse(ComplexF64[1 0; 0 1])
    elseif which == 'X' || which == 1
        return sparse(ComplexF64[0 1; 1 0])
    elseif which == 'Y' || which == 2
        return sparse(ComplexF64[0 -im; im 0])
    elseif which == 'Z' || which == 3
        return sparse(ComplexF64[1 0; 0 -1])
    else
        throw(ArgumentError("Unsupported Pauli label: $which"))
    end
end

end # module Magic

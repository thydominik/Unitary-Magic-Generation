"""
    max_magic (research)

Research script utilities related to searching for high-magic states.

This file was migrated from `src/core/magic/MaxMagic.jl` and refactored so that it:
- Does not execute on import.
- Depends only on core functionality in `src/core/magic/magic.jl`.
- Makes the computational cost explicit.

Warning
The original brute-force coefficient sweep scales as `|C|^(2^n)` where `C` is the coefficient set.
For `n=4` and `|C|=4`, this is `4^16 ≈ 4.3×10^9` candidate states, which is not intended to be
run as-is.
"""
module max_magic

using LinearAlgebra

# Load the core magic routines (non-package project layout).
include(joinpath(@__DIR__, "..", "src", "core", "magic", "magic.jl"))
using .magic

export search_high_magic_states

"""
    search_high_magic_states(; n_qubits=4, coeff_set=[1, im, 1+im, 0], α=2,
                              threshold=-Inf, max_iter=nothing) -> Vector{Vector{ComplexF64}}

Search for states with magic above `threshold` by enumerating coefficient grids.

Parameters
- `n_qubits`: Number of qubits; the state dimension is `2^n_qubits`.
- `coeff_set`: Coefficients used for each computational basis amplitude.
- `α`: Rényi parameter passed to `magic.measure_magic_pure`.
- `threshold`: Keep states whose magic is strictly greater than this value.
- `max_iter`: Optional limit on the number of candidates to test (useful to cap runtime).

Returns
- A vector of normalized state vectors that passed the threshold.

Notes
- This is primarily a *research* utility. For realistic searches, prefer randomized or heuristic
  approaches rather than full enumeration.
"""
function search_high_magic_states(
    ;
    n_qubits::Integer=4,
    coeff_set=(1 + 0im, 0 + 1im, 1 + 1im, 0 + 0im),
    α::Real=2,
    threshold::Real=-Inf,
    max_iter::Union{Nothing, Integer}=nothing,
)
    n_qubits > 0 || throw(ArgumentError("n_qubits must be positive"))
    dim = 2^n_qubits

    strings = magic.generate_all_pauli_strings(n_qubits)
    pauli_ops = magic.pauli_operator_list(strings, n_qubits)

    coeffs = collect(ComplexF64.(coeff_set))
    isempty(coeffs) && throw(ArgumentError("coeff_set must be non-empty"))

    # Enumerate coefficient assignments using base-|coeffs| counting.
    base = length(coeffs)
    total = base^dim
    limit = isnothing(max_iter) ? total : min(total, max_iter)

    out = Vector{Vector{ComplexF64}}()
    state = Vector{ComplexF64}(undef, dim)

    for idx in 0:(limit - 1)
        x = idx
        @inbounds for k in 1:dim
            d = (x % base) + 1
            state[k] = coeffs[d]
            x ÷= base
        end

        nrm = norm(state)
        nrm == 0 && continue
        ψ = state ./ nrm

        m, _ = magic.measure_magic_pure(ψ, pauli_ops, α)
        if m > threshold
            push!(out, copy(ψ))
        end
    end

    return out
end

end # module max_magic

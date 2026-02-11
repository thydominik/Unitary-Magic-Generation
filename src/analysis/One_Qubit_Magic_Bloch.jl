"""
one_qubit_magic_bloch.jl

Compute a 1-qubit magic measure on a Bloch-sphere grid.

This is analysis code (can be slow, produces output files) and is not part of the core library.
It intentionally reuses implementations from src/core.

Notes
- All identifiers and comments are ASCII-only.
- Output paths are relative to this file by default.
"""

using Random

# Optional analysis dependencies.
# using JLD2
# using MAT
# using ProgressBars
# using Plots

include(joinpath(@__DIR__, "..", "core", "magic", "magic.jl"))
using .magic: generate_all_pauli_strings, pauli_operator_list, measure_magic_pure

function compute_magic_bloch(; n_theta::Int=128, n_phi::Int=256, seed::Int=1)
    rng = MersenneTwister(seed)
    Random.seed!(rng)

    n_qubits = 1
    strings = generate_all_pauli_strings(n_qubits)
    pauli_ops = pauli_operator_list(strings, n_qubits)

    theta = collect(range(0.0, pi; length=n_theta))
    phi = collect(range(0.0, 2 * pi; length=n_phi))

    magic_vals = zeros(Float64, length(theta), length(phi))

    # Uncomment ProgressBars usage if installed.
    for i in 1:length(theta)
        for j in 1:length(phi)
            psi = ComplexF64[cos(theta[i] / 2), sin(theta[i] / 2) * exp(im * phi[j])]
            magic_vals[i, j] = measure_magic_pure(psi, pauli_ops, 2)[1]
        end
    end

    return (magic=magic_vals, theta=theta, phi=phi)
end

if abspath(PROGRAM_FILE) == @__FILE__
    out = compute_magic_bloch()
    @info "Computed magic grid with size $(size(out.magic))"

    # Save/plot blocks intentionally omitted by default.
    # Add JLD2/MAT/Plots calls here if you need file output.
end

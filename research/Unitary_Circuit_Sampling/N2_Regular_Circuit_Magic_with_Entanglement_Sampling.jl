"""
N2_Regular_Circuit_Magic_with_Entanglement_Sampling.jl

Same conventions as the N=3 script.
"""

using Random
using ProgressBars
using JLD2
using LinearAlgebra

include(joinpath(@__DIR__, "..", "research_utils.jl"))

current_dir = @__DIR__

include(joinpath(current_dir, "..", "Modules", "Magic.jl"))
include(joinpath(current_dir, "..", "Modules", "Random_Unitaries.jl"))
include(joinpath(current_dir, "..", "Modules", "Entanglement.jl"))

using .Measure_Entanglement
using .Measure_Magic
using .Random_Unitary_Generation

function regular_circuit_magic_with_entanglement_sampling(no_qubits::Int, no_samples::Int, seed::Int)
    Random.seed!(seed)

    psi_0 = 1 / sqrt(2^no_qubits) * ones(ComplexF64, 2^no_qubits)

    pauli_ops = Measure_Magic.PauliOperatorList(Measure_Magic.GenerateAllPauliStrings(no_qubits), no_qubits)

    magic_vals = Vector{Float64}(undef, no_samples)
    entanglement = Vector{Float64}(undef, no_samples)

    for i in ProgressBar(1:no_samples)
        u = Random_Unitary_Generation.Generate_Regular_Unitary_Circuit(no_qubits)
        state = u * psi_0

        rho_left, _ = Measure_Entanglement.reduced_density_matrix(state, [1])
        entanglement[i] = Measure_Entanglement.von_neumann_entropy(rho_left)

        magic_vals[i] = Measure_Magic.MeasureMagic_Pure(state, pauli_ops, 2)[1]
    end

    return magic_vals, entanglement
end

function main()
    save_output = get_bool_env("SAVE_OUTPUT", false)

    no_qubits = 2
    no_samples = 2^20
    seed = 2

    magic_vals, entanglement = regular_circuit_magic_with_entanglement_sampling(no_qubits, no_samples, seed)

    if save_output
        out_jld2 = output_path(
            "RegularUnitaryCircuitMagicEntanglement_N_$(no_qubits)_Samples_$(no_samples)_Seed_$(seed).jld2";
            script_dir=@__DIR__,
            script_file=@__FILE__,
        )
        ensure_parent_dir(out_jld2)
        JLD2.@save out_jld2 magic_vals entanglement no_qubits no_samples seed
        @info "Saved results to $out_jld2"
    else
        @info "SAVE_OUTPUT is disabled; not writing JLD2 output."
    end

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

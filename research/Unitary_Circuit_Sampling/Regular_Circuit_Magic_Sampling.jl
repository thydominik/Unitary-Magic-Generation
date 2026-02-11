"""
Regular_Circuit_Magic_Sampling.jl

Sample magic (and optionally entanglement) from N-qubit random regular unitary circuits.

Important: This is research code. It can take a long time for large N or sample counts.

File output
- Disabled by default.
- Enable with: SAVE_OUTPUT=1 julia Regular_Circuit_Magic_Sampling.jl
- Outputs go to research/output/Unitary_Circuit_Sampling/Regular_Circuit_Magic_Sampling/

Input data
- This script generates data, so it does not require RESEARCH_DATA_DIR.

Notes
- This script still uses legacy include paths ("Modules") from older versions of the repo.
  If you want, we can migrate it to the new src/core API in a follow-up.
"""

using Random
using ProgressBars
using JLD2
using Base.Threads

include(joinpath(@__DIR__, "..", "research_utils.jl"))

current_dir = @__DIR__

# Legacy modules (kept as-is for now).
filepath = joinpath(current_dir, "..", "Modules", "Magic.jl")
include(filepath)
filepath = joinpath(current_dir, "..", "Modules", "Random_Unitaries.jl")
include(filepath)
filepath = joinpath(current_dir, "..", "Modules", "Entanglement.jl")
include(filepath)

using .Measure_Entanglement
using .Measure_Magic
using .Random_Unitary_Generation

function regular_circuit_magic_sampling(no_qubits::Int, no_samples::Int, seed::Int)
    Random.seed!(seed)

    psi_0 = 1 / sqrt(2^no_qubits) * ones(2^no_qubits)

    strings = Measure_Magic.GenerateAllPauliStrings(no_qubits)
    pauli_ops = Measure_Magic.PauliOperatorList(strings, no_qubits)

    magic_vals = Vector{Float64}()
    entanglement = Vector{Vector{Float64}}()
    subsystems = Vector{Vector{Vector{Int}}}()

    for _ in ProgressBar(1:no_samples)
        u = Random_Unitary_Generation.Generate_Regular_Unitary_Circuit(no_qubits)
        state = u * psi_0

        # SVN, subsys = Measure_Entanglement.calculate_entanglement(state, subsystems="Middle")
        # push!(entanglement, SVN)
        # push!(subsystems, subsys)

        push!(magic_vals, Measure_Magic.MeasureMagic_Pure(state, pauli_ops)[1])
    end

    return magic_vals, entanglement, subsystems, psi_0
end

function main()
    save_output = get_bool_env("SAVE_OUTPUT", false)

    n = 7
    seed = 2
    no_samples = 2^20

    magic_vals, entanglement, subsystems, psi_0 = regular_circuit_magic_sampling(n, no_samples, seed)

    if save_output
        out_jld2 = output_path(
            "RegularUnitaryCircuitMagicSampled_N_$(n)_Samples_$(no_samples)_Seed_$(seed).jld2";
            script_dir=@__DIR__,
            script_file=@__FILE__,
        )
        ensure_parent_dir(out_jld2)
        JLD2.@save out_jld2 psi_0 no_samples n seed magic_vals entanglement subsystems
        @info "Saved results to $out_jld2"
    else
        @info "SAVE_OUTPUT is disabled; results are kept in memory only."
    end

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

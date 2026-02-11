"""
N3_Regular_Circuit_Magic_with_Entanglement_Sampling.jl

Sample magic (stabilizer Rényi entropy) and bipartite entanglement cuts for regular random circuits.

File output
- Disabled by default.
- Enable with: SAVE_OUTPUT=1 julia N3_Regular_Circuit_Magic_with_Entanglement_Sampling.jl
- Outputs go to research/output/Unitary_Circuit_Sampling/N3_Regular_Circuit_Magic_with_Entanglement_Sampling/

Notes
- This is research code and can be very slow for large sample counts.
- Uses legacy `research/Modules/*.jl` includes.
"""

using Random
using ProgressBars
using JLD2
using LinearAlgebra

include(joinpath(@__DIR__, "..", "research_utils.jl"))

current_dir = @__DIR__

# Legacy modules
include(joinpath(current_dir, "..", "Modules", "Magic.jl"))
include(joinpath(current_dir, "..", "Modules", "Random_Unitaries.jl"))
include(joinpath(current_dir, "..", "Modules", "Entanglement.jl"))

using .Measure_Entanglement
using .Measure_Magic
using .Random_Unitary_Generation

function regular_circuit_magic_with_entanglement_sampling(no_qubits::Int, no_samples::Int, seed::Int)
    Random.seed!(seed)

    psi_0 = 1 / sqrt(2^no_qubits) * ones(ComplexF64, 2^no_qubits)

    pauli_ops = Dict{Int, Any}()
    for q in 1:no_qubits
        pauli_ops[q] = Measure_Magic.PauliOperatorList(Measure_Magic.GenerateAllPauliStrings(q), q)
    end

    magic_vals = Vector{Float64}(undef, no_samples)
    entanglement = Vector{Vector{Float64}}(undef, no_samples)
    subsystems = Vector{Vector{Vector{Int}}}(undef, no_samples)

    for i in ProgressBar(1:no_samples)
        u = Random_Unitary_Generation.Generate_Regular_Unitary_Circuit(no_qubits)
        state = u * psi_0

        entropies = Float64[]
        subsys = Vector{Int}[]

        for q in 1:(no_qubits - 1)
            left = collect(1:q)
            right = collect((q + 1):no_qubits)

            rho_left, keep_left = Measure_Entanglement.reduced_density_matrix(state, left)
            rho_right, keep_right = Measure_Entanglement.reduced_density_matrix(state, right)

            push!(entropies, Measure_Entanglement.von_neumann_entropy(rho_left))
            push!(subsys, keep_left)

            # keep_right is redundant given keep_left, but we store both for completeness
            push!(subsys, keep_right)

            # (Optional) Mixed-state magic / purity can be added back later if needed.
            _ = rho_right
        end

        entanglement[i] = entropies
        subsystems[i] = subsys
        magic_vals[i] = Measure_Magic.MeasureMagic_Pure(state, pauli_ops[no_qubits], 2)[1]
    end

    return magic_vals, entanglement, subsystems
end

function main()
    save_output = get_bool_env("SAVE_OUTPUT", false)

    no_qubits = 3
    no_samples = 2^20
    seed = 2

    magic_vals, entanglement, subsystems = regular_circuit_magic_with_entanglement_sampling(no_qubits, no_samples, seed)

    if save_output
        out_jld2 = output_path(
            "RegularUnitaryCircuitMagicEntanglement_N_$(no_qubits)_Samples_$(no_samples)_Seed_$(seed).jld2";
            script_dir=@__DIR__,
            script_file=@__FILE__,
        )
        ensure_parent_dir(out_jld2)
        JLD2.@save out_jld2 magic_vals entanglement subsystems no_qubits no_samples seed
        @info "Saved results to $out_jld2"
    else
        @info "SAVE_OUTPUT is disabled; not writing JLD2 output."
    end

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

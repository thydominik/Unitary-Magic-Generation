"""
One_Qubit_Nonstabilising_power.jl

A small 1-qubit experiment: sample random unitaries and measure the magic generated
from each stabiliser input state.

File output
- Disabled by default.
- Enable with: SAVE_OUTPUT=1 julia One_Qubit_Nonstabilising_power.jl
- Outputs go to research/output/Unitary_Circuit_Sampling/One_Qubit_Nonstabilising_power/

Notes
- This script still uses legacy Modules/*.jl includes.
- The plotting section is kept but saving is guarded.
"""

using Random
using ProgressBars
using JLD2
using Base.Threads
using MAT

include(joinpath(@__DIR__, "..", "research_utils.jl"))

current_dir = @__DIR__

# Legacy modules.
filepath = joinpath(current_dir, "..", "Modules", "Magic.jl")
include(filepath)
filepath = joinpath(current_dir, "..", "Modules", "Random_Unitaries.jl")
include(filepath)
filepath = joinpath(current_dir, "..", "Modules", "Entanglement.jl")
include(filepath)

using .Measure_Entanglement
using .Measure_Magic
using .Random_Unitary_Generation

function main()
    save_output = get_bool_env("SAVE_OUTPUT", false)

    seed = 1
    Random.seed!(seed)

    no_samples = 2^22
    no_stab = 6
    no_qubits = 1

    psi = zeros(ComplexF64, 6, 2)
    psi[1, :] = [0, 1]
    psi[2, :] = [1, 0]
    psi[3, :] = [1 / sqrt(2), 1 / sqrt(2)]
    psi[4, :] = [1 / sqrt(2), -1 / sqrt(2)]
    psi[5, :] = [1 / sqrt(2), im / sqrt(2)]
    psi[6, :] = [1 / sqrt(2), -im / sqrt(2)]

    strings = Measure_Magic.GenerateAllPauliStrings(no_qubits)
    pauli_ops = Measure_Magic.PauliOperatorList(strings, no_qubits)

    magic = zeros(Float64, no_samples, no_stab)

    for i in ProgressBar(1:no_samples)
        u = Random_Unitary_Generation.Generate_Regular_Unitary_Circuit(no_qubits)
        for stabstate in 1:no_stab
            psi_0 = psi[stabstate, :]
            state = u * psi_0
            magic[i, stabstate] = Measure_Magic.MeasureMagic_Pure(state, pauli_ops, 2)[1]
        end
    end

    if save_output
        out_jld2 = output_path(
            "UnitaryMagicSampling_AllStabs_N_1_Samples_$(no_samples).jld2";
            script_dir=@__DIR__,
            script_file=@__FILE__,
        )
        ensure_parent_dir(out_jld2)
        JLD2.@save out_jld2 magic
        @info "Saved results to $out_jld2"
    else
        @info "SAVE_OUTPUT is disabled; not writing JLD2 output."
    end

    # Plotting (optional)
    using Plots
    p = plot()
    histogram!(p, magic[:, 1], bins=range(0, log2(3 / 2); length=1000), normalize=:pdf)

    nsp = zeros(Float64, no_samples)
    for i in 1:no_samples
        nsp[i] = sum(magic[i, :]) / no_stab
    end

    histogram!(p, nsp, bins=range(0, log2(3 / 2); length=1000), normalize=:pdf)

    if save_output
        out_pdf = output_path("UnitaryMagicSampling_Histograms.pdf"; script_dir=@__DIR__, script_file=@__FILE__)
        ensure_parent_dir(out_pdf)
        savefig(p, out_pdf)
    end

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

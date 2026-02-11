"""
MagicDistribution_N5.jl

Plot and compare normalised magic distributions for N=5 across circuit depths.

Inputs
- Expects datasets under ENV["RESEARCH_DATA_DIR"].

Outputs
- Saving is disabled by default.
- Enable it with: SAVE_OUTPUT=1 julia MagicDistribution_N5.jl
- Outputs go to research/output/Magic_Distribution_Analysis/MagicDistribution_N5/

Notes
- The original logic is kept, but file paths are now portable.
- This script still assumes the input filenames follow the original naming scheme.
"""

using Statistics
using Plots
using PlotThemes
using LaTeXStrings
using StatsBase
using JLD2

include(joinpath(@__DIR__, "..", "research_utils.jl"))

function main()
    save_output = get_bool_env("SAVE_OUTPUT", false)
    data_dir = get(ENV, "RESEARCH_DATA_DIR", "")
    require_data_dir(data_dir)

    no_qubits = 5
    no_depth = no_qubits + 1
    no_data = no_depth + 1

    # Load data
    data = Dict{Int, Any}()

    reg_path = joinpath(
        data_dir,
        "Random_Unitary_Magic_Generation",
        "RegularUnitaryCircuitMagicSampled_N_$(no_qubits)_Samples_1048576_Seed_1.jld2",
    )
    data[1] = JLD2.load(reg_path)

    for d in 1:no_depth
        bw_path = joinpath(
            data_dir,
            "Random_Unitary_Magic_Generation",
            "BWUnitaryCircuitMagicSampled_N_$(no_qubits)_D_$(d)_Samples_1048576_Seed_1.jld2",
        )
        data[d + 1] = JLD2.load(bw_path)
    end

    # Normalise magic to [0,1]
    m2 = Dict{Int, Vector{Float64}}()
    for i in 1:no_data
        m2[i] = round.(data[i]["Magic"] ./ log((2^no_qubits + 1) / 2); digits=14)
    end

    # Quantiles
    quantiles = Dict{Int, Vector{Float64}}()
    for i in 1:no_data
        quantiles[i] = quantile(m2[i], [0.25, 0.5, 0.75])
    end

    # Histogram
    bin_edges = range(0, 1; length=1001)
    h = Dict{Int, Histogram}()
    for i in 1:no_data
        h[i] = fit(Histogram, m2[i], bin_edges)
        h[i] = StatsBase.normalize(h[i], mode=:pdf)
    end

    counts = Dict{Int, Vector{Float64}}()
    for i in 1:no_data
        counts[i] = h[i].weights
    end

    # Plot: depth vs regular
    p = plot(dpi=400)
    for i in 1:no_data
        if i == 1
            plot!(bin_edges[1:end-1], counts[i], lw=2, label="Regular circuit", legend=:topleft)
        else
            plot!(bin_edges[1:end-1], counts[i], lw=2, label="D = $(i - 1)", legend=:topleft)
        end
    end

    theme(:mute::Symbol;)
    title!(L"Magic Distribution $N = 5$ and $D = 1-6$"; titlefontsize=20)
    xlabel!(L"$\\tilde{M}_2$"; labelfontsize=20)
    ylabel!(L"$\\varrho(\\tilde{M}_2)$"; labelfontsize=20)
    plot!(framestyle=:box)
    plot!(legendfontsize=10)

    if save_output
        out_pdf = output_path("N$(no_qubits)_Depth_vs_reg.pdf"; script_dir=@__DIR__, script_file=@__FILE__)
        out_png = output_path("N$(no_qubits)_Depth_vs_reg.png"; script_dir=@__DIR__, script_file=@__FILE__)
        ensure_parent_dir(out_pdf)
        savefig(p, out_pdf)
        savefig(p, out_png)
    end

    # Plot: regular circuit with quantiles
    p2 = plot(dpi=400)
    plot!(bin_edges[1:end-1], counts[1], lw=2, label="Magic Distribution", legend=:topleft)
    vline!([quantiles[1][1]], color=:red, label="25% = $(round(quantiles[1][1]; digits=5))")
    vline!([quantiles[1][2]], color=:green, label="50% = $(round(quantiles[1][2]; digits=5))")
    vline!([quantiles[1][3]], color=:blue, label="75% = $(round(quantiles[1][3]; digits=5))")
    vline!([mean(m2[1])], linestyle=:dash, label="Average = $(round(mean(m2[1]); digits=4))")
    vline!([median(m2[1])], linestyle=:dash, label="Median = $(round(median(m2[1]); digits=4))")

    theme(:mute::Symbol;)
    title!(L"Magic Distribution $N = 5$"; titlefontsize=20)
    xlabel!(L"$\\tilde{M}_2$"; labelfontsize=20)
    ylabel!(L"$\\varrho(\\tilde{M}_2)$"; labelfontsize=20)
    plot!(framestyle=:box)
    plot!(legendfontsize=10)

    if save_output
        out_pdf = output_path("N$(no_qubits)_RegularCircuit_MagicDistribution.pdf"; script_dir=@__DIR__, script_file=@__FILE__)
        out_png = output_path("N$(no_qubits)_RegularCircuit_MagicDistribution.png"; script_dir=@__DIR__, script_file=@__FILE__)
        ensure_parent_dir(out_pdf)
        savefig(p2, out_pdf)
        savefig(p2, out_png)
    end

    # Sample-size dependence block (kept, but saving guarded)
    samples = zeros(no_data, 20)
    meanv = zeros(no_data, 20)
    medianv = zeros(no_data, 20)
    skewv = zeros(no_data, 20)
    kurtv = zeros(no_data, 20)
    varv = zeros(no_data, 20)

    for keys in 1:no_data, i in 1:20
        samples[keys, i] = 2^i
        meanv[keys, i] = mean(m2[keys][1:Int(samples[keys, i])])
        medianv[keys, i] = median(m2[keys][1:Int(samples[keys, i])])
        skewv[keys, i] = skewness(m2[keys][1:Int(samples[keys, i])])
        kurtv[keys, i] = kurtosis(m2[keys][1:Int(samples[keys, i])])
        varv[keys, i] = var(m2[keys][1:Int(samples[keys, i])])
    end

    if save_output
        out_jld2 = output_path("N$(no_qubits)_Statistics.jld2"; script_dir=@__DIR__, script_file=@__FILE__)
        ensure_parent_dir(out_jld2)
        JLD2.@save out_jld2 counts bin_edges samples meanv medianv skewv kurtv varv
    end

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

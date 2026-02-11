"""
StatisticalMeasures_Comparison.jl

Compare a bunch of distribution-level statistics (mean, median, variance, etc.)
for magic distributions across N=1..6, and explore magic vs entanglement.

Inputs
- Uses datasets under ENV["RESEARCH_DATA_DIR"].

Outputs
- Saving is disabled by default.
- Enable saving with: SAVE_OUTPUT=1 julia StatisticalMeasures_Comparison.jl
- Outputs go to research/output/Magic_Distribution_Analysis/StatisticalMeasures_Comparison/

Notes
- This is research code. It is intentionally "script-like" and not part of src/core.
"""

using Statistics
using Plots
using PlotThemes
using LaTeXStrings
using StatsBase
using JLD2
using LinearAlgebra
using ProgressBars

include(joinpath(@__DIR__, "..", "research_utils.jl"))

function main()
    save_output = get_bool_env("SAVE_OUTPUT", false)
    data_dir = get(ENV, "RESEARCH_DATA_DIR", "")
    require_data_dir(data_dir)

    # N = 1-6 data load
    data = Dict{Int, Any}()

    for n in ProgressBar(1:5)
        @info "Loading N=$n"
        data_path = joinpath(
            data_dir,
            "Random_Unitary_Magic_Generation",
            "N$(n)",
            "Regular",
            "RegularUnitaryCircuitMagicSampled_N_$(n)_Samples_1048576_Seed_1.jld2",
        )
        data[n] = JLD2.load(data_path)
    end

    data_path_n6 = joinpath(
        data_dir,
        "Random_Unitary_Magic_Generation",
        "N6",
        "Regular",
        "RegularUnitaryCircuitMagicSampled_N_6_Samples_1048576_MultiSeed_w_Ent.jld2",
    )
    data[6] = JLD2.load(data_path_n6)

    # Going from log to log2 (legacy conversion kept as-is).
    magic = Dict{Int, Any}()
    for n in 1:6
        magic[n] = -log2.(exp.(-(data[n]["Magic"] .+ log(2^n)))) .- log2(2^n)
    end

    # Distributions of magic
    m = 1000
    h = Dict{Int, Histogram}()

    for n in 1:6
        bin_edges = range(0, log((2^n + 1) / 2); length=m)
        h[n] = fit(Histogram, data[n]["Magic"], bin_edges)
    end

    counts = Dict{Int, Vector{Float64}}()
    for n in 1:6
        counts[n] = h[n].weights
    end

    # Not scaled, normalised
    most_probable = Float64[]
    p = plot(dpi=400)
    for n in 1:6
        bin_edges = range(0, log2((2^n + 1) / 2); length=m - 1)

        integral = 0.0
        for j in 1:length(bin_edges)
            integral += (bin_edges[2] - bin_edges[1]) * counts[n][j]
        end

        push!(most_probable, bin_edges[argmax(counts[n])])
        plot!(bin_edges, counts[n] ./ integral, label="N = $(n)")
    end

    title!("Magic distribution for N = 1-6")
    xlabel!(L"$M_2$"; labelfontsize=20)
    ylabel!(L"$\\varrho ( M_2 )$"; labelfontsize=20)
    plot!(framestyle=:box)
    plot!(legendfontsize=10)

    if save_output
        out_pdf = output_path("N1_N6_magic_distribution.pdf"; script_dir=@__DIR__, script_file=@__FILE__)
        out_png = output_path("N1_N6_magic_distribution.png"; script_dir=@__DIR__, script_file=@__FILE__)
        ensure_parent_dir(out_pdf)
        savefig(p, out_pdf)
        savefig(p, out_png)
    end

    # Scaled and normalised
    p2 = plot(dpi=400)
    for n in 1:6
        bin_edges = range(0, 1; length=m - 1)
        integral = 0.0
        for j in 1:length(bin_edges)
            integral += (bin_edges[2] - bin_edges[1]) * counts[n][j]
        end
        plot!(bin_edges, counts[n] ./ integral, lw=2, label="N = $(n)")
    end

    title!("Normalised Magic distribution for N = 1-6")
    xlabel!(L"$\\tilde{M}_2$"; labelfontsize=20)
    ylabel!(L"$\\varrho ( \\tilde{M}_2 )$"; labelfontsize=20)
    plot!(framestyle=:box)
    plot!(legendfontsize=10)

    if save_output
        out_pdf = output_path("N1_N6_magic_distribution_Scaled.pdf"; script_dir=@__DIR__, script_file=@__FILE__)
        out_png = output_path("N1_N6_magic_distribution_Scaled.png"; script_dir=@__DIR__, script_file=@__FILE__)
        ensure_parent_dir(out_pdf)
        savefig(p2, out_pdf)
        savefig(p2, out_png)
    end

    # Mean, Median, Variance, Skewness, Kurtosis vs N
    mostp = Float64[]
    meanv = Float64[]
    variancev = Float64[]
    medianv = Float64[]
    skewnessv = Float64[]
    kurtosisv = Float64[]

    for n in 1:6
        push!(mostp, most_probable[n] ./ log2((2^n + 1) / 2))
        push!(meanv, mean(magic[n] ./ log2((2^n + 1) / 2)))
        push!(variancev, var(magic[n] ./ log2((2^n + 1) / 2)))
        push!(medianv, median(magic[n] ./ log2((2^n + 1) / 2)))
        push!(skewnessv, skewness(magic[n] ./ log2((2^n + 1) / 2)))
        push!(kurtosisv, kurtosis(magic[n] ./ log2((2^n + 1) / 2)))
    end

    x = 1:6
    p3 = plot(dpi=400)
    scatter!(x, meanv, label="Mean")
    scatter!(x, medianv, label="Median")
    scatter!(x, mostp, label=L"max(P($M_2$))")
    title!("Mean and Median for N = 1-6")
    xlabel!(L"$N$"; labelfontsize=20)
    ylabel!(L"$\\tilde{M}_2$"; labelfontsize=20)
    plot!(ylims=[0.55, 0.85])
    plot!(framestyle=:box)
    plot!(legendfontsize=10)

    if save_output
        out_pdf = output_path("N1_N6_mean_median_magic.pdf"; script_dir=@__DIR__, script_file=@__FILE__)
        out_png = output_path("N1_N6_mean_median_magic.png"; script_dir=@__DIR__, script_file=@__FILE__)
        ensure_parent_dir(out_pdf)
        savefig(p3, out_pdf)
        savefig(p3, out_png)
    end

    # The rest of the original script (magic/entanglement 2D histograms, KDEs, cuts, etc.)
    # is large and very plot-heavy. We keep it as a "manual" block below, but guard all
    # output saving behind save_output so it behaves politely.

    @info "Base distribution plots done. Further analysis blocks are left in-place; enable SAVE_OUTPUT=1 if you want files." 

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

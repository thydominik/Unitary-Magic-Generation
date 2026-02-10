"""
Statistical_Prop_comp.jl

Quick comparison plots of basic statistics (mean/median/skewness/kurtosis)
for the magic distribution across several qubit counts.

This is research code, so it is allowed to have machine-specific data locations.
To make it portable, this script reads the base data directory from:

    ENV["RESEARCH_DATA_DIR"]

Saving figures is off by default. Enable it with:

    SAVE_OUTPUT=1 julia Statistical_Prop_comp.jl

Outputs (when enabled) go to:
- research/output/Magic_Distribution_Analysis/Statistical_Prop_comp/<files>
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

    statdata = Dict{Int, Any}()
    for n in 1:5
        # Note: these are local helper files produced by other scripts.
        statdata[n] = JLD2.load(joinpath(@__DIR__, "N$(n)_Statistics.jld2"))
    end

    data = Dict{Int, Any}()
    for n in 1:5
        @info "Loading N=$n"
        data_path = joinpath(
            data_dir,
            "Random_Unitary_Magic_Generation",
            "RegularUnitaryCircuitMagicSampled_N_$(n)_Samples_1048576_Seed_1.jld2",
        )
        data[n] = JLD2.load(data_path)
    end

    # Normalise magic to [0,1] by dividing by the max possible value log((2^N + 1)/2).
    m2 = Dict{Int, Vector{Float64}}()
    for n in 1:5
        m2[n] = round.(data[n]["Magic"] ./ log((2^n + 1) / 2); digits=14)
    end

    # Distribution comparison
    no_data = 5
    h = Dict{Int, Histogram}()
    bin_edges = range(0, 1; length=1001)
    for n in 1:no_data
        h[n] = fit(Histogram, m2[n], bin_edges)
        h[n] = StatsBase.normalize(h[n], mode=:pdf)
    end

    counts = Dict{Int, Vector{Float64}}()
    for n in 1:no_data
        counts[n] = h[n].weights
    end

    p = plot(dpi=400)
    for n in 1:5
        plot!(bin_edges[1:end-1], counts[n], lw=2, label="N = $(n)", legend=:topleft)
    end

    theme(:mute::Symbol;)
    title!(L"Magic Distribution $N = 1-5$"; titlefontsize=20)
    xlabel!(L"$\\tilde{M}_2$"; labelfontsize=20)
    ylabel!(L"$\\varrho(\\tilde{M}_2)$"; labelfontsize=20)
    plot!(framestyle=:box)
    plot!(legendfontsize=10)
    vline!([mean(m2[1]), mean(m2[2]), mean(m2[3]), mean(m2[4]), mean(m2[5])], color="black", label="averages")

    if save_output
        out_pdf = output_path("n_1_5_normalised.pdf"; script_dir=@__DIR__, script_file=@__FILE__)
        out_png = output_path("n_1_5_normalised.png"; script_dir=@__DIR__, script_file=@__FILE__)
        ensure_parent_dir(out_pdf)
        savefig(p, out_pdf)
        savefig(p, out_png)
    else
        @info "SAVE_OUTPUT is disabled; not writing figures."
    end

    # Mean and median
    mean_magic = Float64[]
    for n in 1:5
        if n == 1 || n == 2
            push!(mean_magic, statdata[n]["Mean"][end])
        else
            push!(mean_magic, maximum(statdata[n]["Mean"]))
        end
    end

    nvals = [1, 2, 3, 4, 5]
    p2 = plot(dpi=400)
    scatter!(nvals, [mean(m2[1]), mean(m2[2]), mean(m2[3]), mean(m2[4]), mean(m2[5])], label="Mean")
    scatter!(nvals, [median(m2[1]), median(m2[2]), median(m2[3]), median(m2[4]), median(m2[5])], label="Median")
    plot!(framestyle=:box)
    plot!(legendfontsize=10)
    xlabel!(L"N"; labelfontsize=20)
    ylabel!("Statistical Vars"; labelfontsize=20)

    if save_output
        out_pdf = output_path("n_1_5_mean_median.pdf"; script_dir=@__DIR__, script_file=@__FILE__)
        out_png = output_path("n_1_5_mean_median.png"; script_dir=@__DIR__, script_file=@__FILE__)
        ensure_parent_dir(out_pdf)
        savefig(p2, out_pdf)
        savefig(p2, out_png)
    end

    p3 = plot(dpi=400)
    scatter!(nvals, [skewness(m2[1]), skewness(m2[2]), skewness(m2[3]), skewness(m2[4]), skewness(m2[5])], label="Skewness")
    scatter!(nvals, [kurtosis(m2[1]), kurtosis(m2[2]), kurtosis(m2[3]), kurtosis(m2[4]), kurtosis(m2[5])], label="Kurtosis")
    plot!(framestyle=:box)
    plot!(legendfontsize=10)
    xlabel!(L"N"; labelfontsize=20)
    ylabel!("Statistical Vars"; labelfontsize=20)

    if save_output
        out_pdf = output_path("n_1_5_skew_kurt.pdf"; script_dir=@__DIR__, script_file=@__FILE__)
        out_png = output_path("n_1_5_skew_kurt.png"; script_dir=@__DIR__, script_file=@__FILE__)
        ensure_parent_dir(out_pdf)
        savefig(p3, out_pdf)
        savefig(p3, out_png)
    end

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

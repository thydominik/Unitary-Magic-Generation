"""
MagicDistribution_N2.jl

See MagicDistribution_N4.jl for conventions.
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

    no_qubits = 2
    no_depth = no_qubits + 1
    no_data = no_depth + 1

    data = Dict{Int, Any}()

    reg_path = joinpath(data_dir, "Random_Unitary_Magic_Generation",
                        "RegularUnitaryCircuitMagicSampled_N_$(no_qubits)_Samples_1048576_Seed_1.jld2")
    data[1] = JLD2.load(reg_path)

    for d in 1:no_depth
        bw_path = joinpath(data_dir, "Random_Unitary_Magic_Generation",
                           "BWUnitaryCircuitMagicSampled_N_$(no_qubits)_D_$(d)_Samples_1048576_Seed_1.jld2")
        data[d + 1] = JLD2.load(bw_path)
    end

    m2 = Dict{Int, Vector{Float64}}()
    for i in 1:no_data
        m2[i] = round.(data[i]["Magic"] ./ log((2^no_qubits + 1) / 2); digits=14)
    end

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

    p = plot(dpi=400)
    for i in 1:no_data
        label = (i == 1) ? "Regular circuit" : "D = $(i - 1)"
        plot!(bin_edges[1:end-1], counts[i], lw=2, label=label, legend=:topleft)
    end

    theme(:mute::Symbol;)
    title!(L"Magic Distribution $N = 2$ and $D = 1-3$"; titlefontsize=20)
    xlabel!(L"$\\tilde{M}_2$"; labelfontsize=20)
    ylabel!(L"$\\varrho(\\tilde{M}_2)$"; labelfontsize=20)
    plot!(framestyle=:box)
    plot!(legendfontsize=10)

    if save_output
        out_pdf = output_path("N2_Depth_vs_reg.pdf"; script_dir=@__DIR__, script_file=@__FILE__)
        out_png = output_path("N2_Depth_vs_reg.png"; script_dir=@__DIR__, script_file=@__FILE__)
        ensure_parent_dir(out_pdf)
        savefig(p, out_pdf)
        savefig(p, out_png)
    end

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

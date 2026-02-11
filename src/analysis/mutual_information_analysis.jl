"""
mutual_information_analysis.jl

Example script to compute and visualize mutual information for stored datasets.

This is analysis code (may use plotting and file I/O) and is not part of the core library.
It intentionally reuses implementations from src/core.

Notes
- All identifiers and comments are ASCII-only.
- Update `data_paths` below to match your machine.
"""

using Random

# Optional dependencies for analysis. Install as needed.
# using JLD2
# using Plots
# using ProgressBars

include(joinpath(@__DIR__, "..", "core", "mutual_information", "mutual_information.jl"))
using .mutual_information: mutual_information_histogram

"""Build a data filename from a base directory and a naming pattern."""
function build_data_filename(base_dir::AbstractString, pattern::AbstractString, n_qubits::Int, n_samples::Int)::String
    name = replace(replace(pattern, "{N}" => string(n_qubits)), "{SAMPLES}" => string(n_samples))
    return joinpath(base_dir, name)
end

"""Compute MI vs sample size for exponentially increasing sample sizes."""
function mi_vs_sample_size(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; n_trials::Int=14, n_bins::Int=64, seed::Int=1234)
    length(x) == length(y) || throw(DimensionMismatch("x and y must have the same length"))

    rng = MersenneTwister(seed)
    mi_vals = Float64[]
    sample_sizes = Int[]

    for m in 0:(n_trials - 1)
        s = 2^m
        s > length(x) && break

        idx = randperm(rng, length(x))[1:s]
        mi = mutual_information_histogram(x[idx], y[idx]; n_bins=n_bins)
        push!(mi_vals, mi)
        push!(sample_sizes, s)
    end

    return mi_vals, sample_sizes
end

if abspath(PROGRAM_FILE) == @__FILE__
    # ---------------------------------------------------------------------
    # Configure your data here.
    # ---------------------------------------------------------------------
    data_paths = Dict(
        "base_dir" => "/path/to/data",
        "file_pattern" => "RegularUnitaryCircuitMagicSampled_N_{N}_Samples_{SAMPLES}.jld2",
        "n_qubits" => 2:10,
        "sample_counts" => [
            1048576,
            1048576,
            1048576,
            1048576,
            1048576,
            1048576,
            442368,
            1048576,
            196608,
        ],
    )

    @info "This script is a template. Uncomment JLD2/Plots imports and implement loading for your dataset."
    @info "Core MI estimator available as mutual_information.mutual_information_histogram"
end

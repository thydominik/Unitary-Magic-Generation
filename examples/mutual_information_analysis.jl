"""
    mutual_information_analysis.jl

Example script demonstrating the refactored mutual information analysis pipeline.

This script shows how to:
1. Load data from JLD2 files
2. Compute mutual information with different sample sizes
3. Analyze the sample size dependence
4. Generate visualization plots

# Usage
```julia
julia> include("examples/mutual_information_analysis.jl")
```
"""

using Plots
using JLD2
using ProgressBars
using Random
using StatsBase
include("../src/UnitaryMagic.jl")
using .UnitaryMagic: MI, MIn

"""
    load_circuit_data(n_qubits::Int, samples::Vector{Int})::Dict

Load unitary circuit data from JLD2 files.

# Arguments
- `n_qubits::Int`: Range of qubit numbers to load
- `samples::Vector{Int}`: Sample counts for each qubit configuration

# Returns
- `Dict`: Dictionary mapping qubit count to loaded data
"""
function load_circuit_data(n_qubits::UnitRange, samples::Vector{Int})::Dict
    data = Dict()
    for (idx, i) in enumerate(n_qubits)
        filename = "/Volumes/SSD_Szmbthy/PhD_Data/Unitary_Circuits_Magic_and_Entanglement/Regular/RegularUnitaryCircuitMagicSampled_N_$(i)_Samples_$(samples[idx]).jld2"
        if isfile(filename)
            data[i] = load(filename)
            @info "Loaded data for N=$i with $(samples[idx]) samples"
        else
            @warn "File not found: $filename"
        end
    end
    return data
end

"""
    analyze_sample_size_dependence(x::Vector, y::Vector, n_trials::Int=14)::Tuple

Analyze how mutual information depends on sample size.

# Arguments
- `x::Vector`: First data array (e.g., magic)
- `y::Vector`: Second data array (e.g., entanglement)
- `n_trials::Int`: Number of exponentially spaced samples (default: 14)

# Returns
- Tuple of (MI_values::Vector, sample_sizes::Vector)

# Notes
- Sample sizes scale as 2^m for m = 0 to n_trials-1
- Uses fixed bin count of 2^12 for MI computation
"""
function analyze_sample_size_dependence(
    x::Vector, 
    y::Vector, 
    n_trials::Int = 14,
    n_bins::Int = 2^12
)::Tuple{Vector{Float64}, Vector{Int}}
    
    MI_values = Float64[]
    sample_sizes = Int[]
    
    for m in ProgressBar(0:n_trials-1)
        sample_size = 2^m
        if sample_size <= length(x)
            selected_indices = randperm(length(x))[1:sample_size]
            mutual_info, _ = MI(x[selected_indices], y[selected_indices], n_bins)
            push!(MI_values, mutual_info)
            push!(sample_sizes, sample_size)
        end
    end
    
    return MI_values, sample_sizes
end

"""
    plot_sample_size_scaling(results::Dict, n_qubits::Vector{Int}, output_file::String)

Plot mutual information as function of sample size for all qubit counts.
"""
function plot_sample_size_scaling(
    results::Dict, 
    n_qubits::Vector{Int}, 
    output_file::String = "MI_sample_size_scaling.png"
)
    p = plot(
        xlabel = "Sample Size - log₂(M)",
        ylabel = "Mutual Information",
        legend = true,
        box = :on
    )
    
    # Color gradient from blue to red
    n_colors = length(results)
    colors = [RGB(0, 0, 1) * (1 - (i-1)/(n_colors-1)) + RGB(1, 0, 0) * ((i-1)/(n_colors-1)) 
              for i in 1:n_colors]
    
    for (idx, n) in enumerate(n_qubits)
        if haskey(results, n)
            mi_vals, sample_sizes = results[n]
            scatter!(p, log2.(sample_sizes), mi_vals, 
                    label = "N = $n", color = colors[idx], markersize = 5)
        end
    end
    
    savefig(p, output_file)
    @info "Plot saved to $output_file"
    return p
end

"""
    plot_mi_vs_qubit_count(results::Dict, sample_size_exp::Int)

Plot mutual information as function of qubit count for fixed sample size.
"""
function plot_mi_vs_qubit_count(
    results::Dict, 
    sample_size_exp::Int = 20,  # 2^20 samples
    output_file::String = "MI_qubit_dependence.png"
)
    p = plot(
        xlabel = "Number of Qubits (N)",
        ylabel = "Mutual Information",
        legend = true,
        box = :on
    )
    
    n_qubits = sort(collect(keys(results)))
    mi_at_target = Float64[]
    
    for n in n_qubits
        mi_vals, sample_sizes = results[n]
        # Find MI at sample size closest to 2^sample_size_exp
        target_size = 2^sample_size_exp
        if any(abs.(sample_sizes .- target_size) .< 1)
            idx = argmin(abs.(sample_sizes .- target_size))
            push!(mi_at_target, mi_vals[idx])
        end
    end
    
    scatter!(p, n_qubits[1:length(mi_at_target)], mi_at_target, 
            color = :blue, markersize = 8, label = "MI at MS = 2^$sample_size_exp")
    
    savefig(p, output_file)
    @info "Plot saved to $output_file"
    return p
end


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    @info "Starting mutual information analysis pipeline..."
    
    # Configuration
    N_QUBITS = 2:10
    SAMPLE_COUNTS = [1048576, 1048576, 1048576, 1048576, 1048576, 1048576, 442368, 1048576, 196608]
    
    # Load data
    @info "Loading circuit data..."
    data = load_circuit_data(N_QUBITS, SAMPLE_COUNTS)
    
    if !isempty(data)
        # Analyze sample size dependence for each qubit count
        @info "Computing mutual information for different sample sizes..."
        results = Dict()
        
        for n in ProgressBar(sort(collect(keys(data))))
            x = data[n]["Magic"]
            y = data[n]["Svn"]
            results[n] = analyze_sample_size_dependence(x, y)
        end
        
        # Save results
        @save "mi_analysis_results.jld2" results
        @info "Results saved to mi_analysis_results.jld2"
        
        # Generate plots
        plot_sample_size_scaling(results, sort(collect(keys(results))))
        plot_mi_vs_qubit_count(results, 20)
    else
        @warn "No data files found. Please adjust file paths in load_circuit_data()."
    end
end

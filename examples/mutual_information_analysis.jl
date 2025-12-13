"""
    mutual_information_analysis.jl

Example script demonstrating the complete MI analysis pipeline.

This script shows best practices for:
- Importing functions from UnitaryMagic module (NOT duplicating code)
- Separating data configuration from logic
- Computing mutual information with various sample sizes
- Analyzing sample size effects on MI estimation
- Generating publication-quality visualizations

# Usage

```julia
julia> include("examples/mutual_information_analysis.jl")
```

# Configuration

Update the DATA_PATHS dictionary below to point to your actual data files.

# Output

Generates:
- `mi_analysis_results.jld2`: Cached MI computation results
- `MI_sample_size_scaling.png`: Plot of MI vs sample size
- `MI_qubit_dependence.png`: Plot of MI vs qubit count
"""

# ============================================================================
# MODULE IMPORTS - ALL FUNCTIONS CALLED FROM UnitaryMagic MODULE
# ============================================================================

using Plots        # For visualization
using JLD2         # For data I/O
using ProgressBars # For progress tracking
using Random       # For random sampling
using StatsBase    # For statistics

# Import the main package module
include("../src/UnitaryMagic.jl")
using .UnitaryMagic: MI, MIn  # Import FUNCTIONS FROM module, don't redefine them!

# ============================================================================
# CONFIGURATION - SEPARATE FROM LOGIC
# ============================================================================

"""
DATA_PATHS

Configuration dictionary for data loading paths and parameters.

This separates data configuration from analysis logic, making it easy to:
- Switch between different datasets
- Run on different machines
- Test with synthetic data
"""
const DATA_PATHS = Dict(
    # Data directory (modify this for your system)
    "base_dir" => "/Volumes/SSD_Szmbthy/PhD_Data/Unitary_Circuits_Magic_and_Entanglement/Regular/",
    
    # File naming pattern
    "file_pattern" => "RegularUnitaryCircuitMagicSampled_N_{N}_Samples_{SAMPLES}.jld2",
    
    # Qubit range to analyze
    "n_qubits" => 2:10,
    
    # Sample counts for each qubit configuration
    "sample_counts" => [
        1048576,  # N=2:  2^20 samples
        1048576,  # N=3:  2^20 samples
        1048576,  # N=4:  2^20 samples
        1048576,  # N=5:  2^20 samples
        1048576,  # N=6:  2^20 samples
        1048576,  # N=7:  2^20 samples
        442368,   # N=8:  ~2^19 samples
        1048576,  # N=9:  2^20 samples
        196608,   # N=10: ~2^18 samples
    ]
)

"""
MI_PARAMS

Configuration parameters for mutual information computation.
"""
const MI_PARAMS = Dict(
    # Number of histogram bins (controls resolution)
    "n_bins_default" => 2^12,     # 4096 bins - good balance
    
    # Number of sample size trials (samples as 2^m for m=0..n_trials-1)
    "n_trials" => 14,             # Tests 2^0, 2^1, ..., 2^13 samples
    
    # Target sample size for qubit dependence plot (2^20 = ~1M samples)
    "target_sample_size_exp" => 20,
)

"""
OUTPUT_PATHS

Configuration for output files and visualizations.
"""
const OUTPUT_PATHS = Dict(
    "results_file" => "mi_analysis_results.jld2",
    "plot_sample_size" => "MI_sample_size_scaling.png",
    "plot_qubit_dep" => "MI_qubit_dependence.png",
)

# ============================================================================
# PRIVATE HELPER FUNCTIONS
# ============================================================================

"""
    _build_data_filename(n_qubits::Int, n_samples::Int)::String

Internal: Build full path to data file based on configuration.

# Arguments
- `n_qubits::Int`: Number of qubits (N parameter)
- `n_samples::Int`: Number of samples (SAMPLES parameter)

# Returns
- `String`: Full path to data file
"""
function _build_data_filename(n_qubits::Int, n_samples::Int)::String
    base = DATA_PATHS["base_dir"]
    pattern = DATA_PATHS["file_pattern"]
    
    # Replace placeholders in pattern
    filename = replace(
        replace(pattern, "{N}" => n_qubits),
        "{SAMPLES}" => n_samples
    )
    
    return joinpath(base, filename)
end

"""
    _validate_data_dimensions(x::AbstractVector, y::AbstractVector)::Bool

Internal: Validate that data arrays have compatible dimensions.
"""
function _validate_data_dimensions(x::AbstractVector, y::AbstractVector)::Bool
    if length(x) != length(y)
        @warn "Data length mismatch: x=$(length(x)), y=$(length(y))"
        return false
    end
    return true
end

# ============================================================================
# PUBLIC FUNCTIONS
# ============================================================================

"""
    load_circuit_data(
        n_qubits_range::UnitRange, 
        sample_counts::Vector{Int}
    )::Dict

Load unitary circuit data from JLD2 files.

This function loads pre-computed quantum circuit data for MI analysis.
Data should be in JLD2 format with "Magic" and "Svn" keys containing
magic content and entanglement measures.

# Arguments

- `n_qubits_range::UnitRange`: Range of qubit numbers to load (e.g., 2:10)
  - Corresponds to circuit size
  - Used in file path construction

- `sample_counts::Vector{Int}`: Sample counts for each qubit configuration
  - Length must match range length
  - sample_counts[i] is used for n_qubits_range[i]
  - Enables varying sample sizes across configurations

# Returns

- `Dict{Int, Dict}`: Dictionary mapping qubit count to loaded data
  - Keys: qubit numbers (2, 3, 4, ...)
  - Values: Dict with keys "Magic" and "Svn" (at minimum)
  - Empty if no files found (check @warn messages)

# Data Format

Expected JLD2 file structure:
```julia
jld2_file = {
    "Magic" => [m₁, m₂, ..., m_N],       # Magic content of each sample
    "Svn" => [s₁, s₂, ..., s_N],         # Von Neumann entropy of each sample
    "CustomData" => [...],                # Any other fields
}
```

# Example

```julia
data = load_circuit_data(2:4, [1000000, 500000, 250000])
# Loads:
# - N=2 with 1M samples
# - N=3 with 500k samples  
# - N=4 with 250k samples

magic_N3 = data[3]["Magic"]  # Access magic data for N=3
svn_N3 = data[3]["Svn"]      # Access entanglement data for N=3
```

# Implementation Notes

- Non-blocking: Silently skips missing files with @warn message
- Logs successful loads with @info messages
- File path built using DATA_PATHS configuration
- No data validation performed (caller should validate)

# Performance

- Time: O(N_files × file_size) - mostly I/O limited
- Memory: O(total_samples) - all data loaded into RAM
- For large datasets, may need memory management
"""
function load_circuit_data(
    n_qubits_range::UnitRange, 
    sample_counts::Vector{Int}
)::Dict{Int, Dict}
    
    # ========================================================================
    # INPUT VALIDATION
    # ========================================================================
    if length(n_qubits_range) != length(sample_counts)
        error(
            "Dimension mismatch: $(length(n_qubits_range)) qubit configs " *
            "but $(length(sample_counts)) sample counts provided"
        )
    end
    
    # ========================================================================
    # LOAD DATA
    # ========================================================================
    data = Dict{Int, Dict}()
    
    for (idx, n_qubit) in enumerate(n_qubits_range)
        # Build full file path
        n_sample = sample_counts[idx]
        filepath = _build_data_filename(n_qubit, n_sample)
        
        # Try to load file
        if isfile(filepath)
            try
                loaded_data = load(filepath)
                data[n_qubit] = loaded_data
                @info "Successfully loaded N=$n_qubit with $n_sample samples"
            catch err
                @warn "Error loading $filepath: $err"
            end
        else
            @warn "File not found: $filepath"
        end
    end
    
    return data
end

"""
    analyze_sample_size_dependence(
        x::AbstractVector, 
        y::AbstractVector, 
        n_trials::Int=14,
        n_bins::Int=2^12
    )::Tuple{Vector{Float64}, Vector{Int}}

Analyze how mutual information depends on sample size.

This function computes MI for exponentially-increasing sample sizes to study
how MI estimation converges as more samples are available.

# Arguments

- `x::AbstractVector`: First data array (e.g., magic content)
  - Length should be large (1000+)
  - Will be subsampled for different sample sizes

- `y::AbstractVector`: Second data array (e.g., entanglement)
  - Must have same length as x
  - Corresponding samples from x

- `n_trials::Int`: Number of sample sizes to test (default: 14)
  - Tests 2^0, 2^1, 2^2, ..., 2^(n_trials-1) samples
  - More trials = finer resolution but slower

- `n_bins::Int`: Number of histogram bins for MI computation (default: 2^12)
  - See MI() documentation for guidance
  - Fixed across all sample sizes for consistency

# Returns

- Tuple of (mi_values, sample_sizes):
  - `mi_values::Vector{Float64}`: MI at each sample size
  - `sample_sizes::Vector{Int}`: Corresponding sample counts (2^0, 2^1, ...)
  - Both vectors have same length

# Example

```julia
x = randn(100000)
y = 0.7 .* x .+ 0.3 .* randn(100000)

mi_vals, sample_sizes = analyze_sample_size_dependence(x, y, n_trials=12)

# Results
println("Sample sizes: \$sample_sizes")
println("MI values: \$mi_vals")
println("Convergence: MI grows from $(mi_vals[1]) to $(mi_vals[end])")
```

# Algorithm

For m = 0 to n_trials-1:
1. Sample size = 2^m
2. Randomly select 2^m samples from data
3. Compute MI using MI() function from UnitaryMagic module
4. Store result in arrays
5. Display progress with ProgressBar

# Performance Notes

- Time: O(n_trials × n_bins²) for MI computations
  - Plus O(n_trials × 2^n_trials) for random sampling
- Memory: O(2^n_trials) for storing subsampled data
- Each MI computation takes ~0.1-1.0 seconds depending on n_bins

# Convergence Analysis

MI estimates should stabilize as sample size increases:
- Few samples: MI varies widely
- Many samples: MI converges to true value
- Rate of convergence depends on:
  - Distribution complexity
  - Dimensionality
  - Dependence strength

"""
function analyze_sample_size_dependence(
    x::AbstractVector, 
    y::AbstractVector, 
    n_trials::Int = MI_PARAMS["n_trials"],
    n_bins::Int = MI_PARAMS["n_bins_default"]
)::Tuple{Vector{Float64}, Vector{Int}}
    
    # ========================================================================
    # INPUT VALIDATION
    # ========================================================================
    if !_validate_data_dimensions(x, y)
        error("Invalid data dimensions")
    end
    
    if n_trials < 1
        error("n_trials must be >= 1, got $n_trials")
    end
    
    # ========================================================================
    # ANALYZE SAMPLE SIZE DEPENDENCE
    # ========================================================================
    mi_values = Float64[]
    sample_sizes = Int[]
    
    # Loop through exponential sample sizes
    for m in ProgressBar(0:n_trials-1)
        sample_size = 2^m
        
        # Skip if sample size exceeds available data
        if sample_size > length(x)
            @debug "Skipping m=$m (2^$m > $(length(x)) samples available)"
            break
        end
        
        # Randomly select subset of data
        selected_indices = randperm(length(x))[1:sample_size]
        x_subset = x[selected_indices]
        y_subset = y[selected_indices]
        
        # IMPORTANT: Call MI() function FROM UnitaryMagic module
        # NOT a local reimplementation
        mi_result, _ = MI(x_subset, y_subset, n_bins)
        
        # Store results
        push!(mi_values, mi_result)
        push!(sample_sizes, sample_size)
    end
    
    return mi_values, sample_sizes
end

"""
    plot_sample_size_scaling(
        results::Dict{Int, Tuple}, 
        n_qubits::Vector{Int},
        output_file::String="MI_sample_size_scaling.png"
    )::Plots.Plot

Plot mutual information as function of sample size for all qubit counts.

Creates a multi-curve plot showing how MI estimation converges with sample
size for different numbers of qubits. Useful for publication figures.

# Arguments

- `results::Dict{Int, Tuple}`: Results from analyze_sample_size_dependence()
  - Keys: qubit numbers
  - Values: (mi_values, sample_sizes) tuples

- `n_qubits::Vector{Int}`: Qubit numbers to plot (in order)
  - Determines colors and labels
  - Should match keys in results dict

- `output_file::String`: Path for saved PNG (default: "MI_sample_size_scaling.png")
  - Can include directory path
  - PNG format suitable for papers

# Returns

- `Plots.Plot`: The generated plot object
  - Can be further customized if needed
  - Also saved to output_file

# Example

```julia
results = Dict(
    2 => (mi_vals_2, sample_sizes_2),
    3 => (mi_vals_3, sample_sizes_3),
    4 => (mi_vals_4, sample_sizes_4),
)

plot_sample_size_scaling(results, [2, 3, 4], "output/scaling.png")
```

# Visualization Features

- X-axis: log₂ scale of sample size (0 to n_trials-1)
- Y-axis: Mutual information in bits
- Different colors for each qubit count (blue→red gradient)
- Legend with qubit labels
- Grid and box for readability
"""
function plot_sample_size_scaling(
    results::Dict{Int, Tuple}, 
    n_qubits::Vector{Int},
    output_file::String = OUTPUT_PATHS["plot_sample_size"]
)::Plots.Plot
    
    # ========================================================================
    # CREATE PLOT FIGURE
    # ========================================================================
    p = plot(
        xlabel = "Sample Size - log\u2082(M)",
        ylabel = "Mutual Information (bits)",
        title = "MI Convergence vs Sample Size",
        legend = :bottomright,
        box = :on,
        size = (800, 600)
    )
    
    # ========================================================================
    # ADD DATA CURVES WITH COLOR GRADIENT
    # ========================================================================
    # Create color gradient from blue (low N) to red (high N)
    n_colors = length(results)
    colors = if n_colors > 1
        [RGB(0, 0, 1) * (1 - (i-1)/(n_colors-1)) + RGB(1, 0, 0) * ((i-1)/(n_colors-1)) 
         for i in 1:n_colors]
    else
        [RGB(0, 0, 1)]
    end
    
    # Plot each qubit configuration
    for (idx, n) in enumerate(n_qubits)
        if haskey(results, n)
            mi_vals, sample_sizes = results[n]
            sample_size_log2 = log2.(sample_sizes)
            
            scatter!(
                p, 
                sample_size_log2, 
                mi_vals,
                label = "N = $n",
                color = colors[idx],
                markersize = 6,
                markerstrokewidth = 1,
                alpha = 0.8
            )
        else
            @warn "No results for N=$n"
        end
    end
    
    # ========================================================================
    # SAVE AND RETURN
    # ========================================================================
    savefig(p, output_file)
    @info "Plot saved to $output_file"
    return p
end

"""
    plot_mi_vs_qubit_count(
        results::Dict{Int, Tuple}, 
        sample_size_exp::Int=20,
        output_file::String="MI_qubit_dependence.png"
    )::Plots.Plot

Plot mutual information as function of qubit count for fixed sample size.

Extracts MI values at a specific sample size across different qubit counts,
showing system-size dependence of measured dependence.

# Arguments

- `results::Dict{Int, Tuple}`: Results from analyze_sample_size_dependence()
- `sample_size_exp::Int`: Target sample size exponent (default: 20 means 2^20 ≈ 1M)
- `output_file::String`: Output PNG file path

# Returns

- `Plots.Plot`: The generated plot

# Example

```julia
plot_mi_vs_qubit_count(results, 20)  # MI at 2^20 = 1M samples
```
"""
function plot_mi_vs_qubit_count(
    results::Dict{Int, Tuple}, 
    sample_size_exp::Int = MI_PARAMS["target_sample_size_exp"],
    output_file::String = OUTPUT_PATHS["plot_qubit_dep"]
)::Plots.Plot
    
    # ========================================================================
    # EXTRACT TARGET SAMPLE SIZE DATA
    # ========================================================================
    n_qubits = sort(collect(keys(results)))
    mi_at_target = Float64[]
    n_valid = Int[]
    
    target_size = 2^sample_size_exp
    
    for n in n_qubits
        mi_vals, sample_sizes = results[n]
        
        # Find MI value at sample size closest to target
        if !isempty(sample_sizes)
            idx = argmin(abs.(sample_sizes .- target_size))
            push!(mi_at_target, mi_vals[idx])
            push!(n_valid, n)
        end
    end
    
    # ========================================================================
    # CREATE PLOT
    # ========================================================================
    p = plot(
        xlabel = "Number of Qubits (N)",
        ylabel = "Mutual Information (bits)",
        title = "MI vs Qubit Count at 2^$(sample_size_exp) Samples",
        legend = :topleft,
        box = :on,
        size = (800, 600)
    )
    
    # Plot data points and trend
    scatter!(
        p,
        n_valid,
        mi_at_target,
        color = :blue,
        markersize = 8,
        markerstrokewidth = 1,
        label = "MI at 2^$(sample_size_exp) samples"
    )
    
    # Try to add trend line if enough points
    if length(n_valid) >= 3
        try
            # Linear fit through points
            A = [ones(length(n_valid)) n_valid]
            coeffs = A \ mi_at_target
            trend_line = A * coeffs
            plot!(
                p,
                n_valid,
                trend_line,
                color = :red,
                linestyle = :dash,
                label = "Linear trend",
                alpha = 0.7
            )
        catch err
            @debug "Could not fit trend line: $err"
        end
    end
    
    # ========================================================================
    # SAVE AND RETURN
    # ========================================================================
    savefig(p, output_file)
    @info "Plot saved to $output_file"
    return p
end

# ============================================================================
# MAIN EXECUTION
# ============================================================================

"""
Main execution block - only runs if this script is directly executed.

Does NOT run if this file is included from another script.
"""
if abspath(PROGRAM_FILE) == @__FILE__
    @info "="^80
    @info "Starting Mutual Information Analysis Pipeline"
    @info "="^80
    
    # ========================================================================
    # STEP 1: LOAD DATA
    # ========================================================================
    @info "\nStep 1: Loading circuit data..."
    data = load_circuit_data(
        DATA_PATHS["n_qubits"],
        DATA_PATHS["sample_counts"]
    )
    
    if isempty(data)
        @error "No data loaded. Check DATA_PATHS configuration."
        exit(1)
    end
    
    @info "Loaded data for $(length(data)) qubit configurations"
    
    # ========================================================================
    # STEP 2: ANALYZE SAMPLE SIZE DEPENDENCE
    # ========================================================================
    @info "\nStep 2: Computing MI for different sample sizes..."
    results = Dict{Int, Tuple}()
    
    for n in ProgressBar(sort(collect(keys(data))))
        @info "  Processing N=$n..."
        x = data[n]["Magic"]
        y = data[n]["Svn"]
        
        # Validate data
        if !_validate_data_dimensions(x, y)
            @warn "Skipping N=$n due to data dimension mismatch"
            continue
        end
        
        # IMPORTANT: Call analyze_sample_size_dependence which uses MI() FROM module
        results[n] = analyze_sample_size_dependence(
            x, y,
            MI_PARAMS["n_trials"],
            MI_PARAMS["n_bins_default"]
        )
    end
    
    @info "Computed MI for $(length(results)) configurations"
    
    # ========================================================================
    # STEP 3: SAVE RESULTS
    # ========================================================================
    @info "\nStep 3: Saving analysis results..."
    output_file = OUTPUT_PATHS["results_file"]
    @save output_file results
    @info "Results saved to $output_file"
    
    # ========================================================================
    # STEP 4: GENERATE VISUALIZATIONS
    # ========================================================================
    @info "\nStep 4: Generating visualizations..."
    plot_sample_size_scaling(
        results,
        sort(collect(keys(results)))
    )
    
    plot_mi_vs_qubit_count(
        results,
        MI_PARAMS["target_sample_size_exp"]
    )
    
    @info "\n" * "="^80
    @info "Analysis complete!"
    @info "="^80
end

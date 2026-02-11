"""
Gaussian mutual information experiments.

This script was moved out of src/core because it performs plotting and long-running experiments.
All identifiers and comments are ASCII-only.
"""

using LinearAlgebra
using Random

# Local includes (project-style repository)
include(joinpath(@__DIR__, "..", "..", "core", "utilities", "utilities.jl"))
include(joinpath(@__DIR__, "..", "..", "core", "mutual_information", "mutual_information.jl"))

using .utilities
using .mutual_information

# Example: compare numerical MI for correlated Gaussian data against analytic value.
function gaussian_mi_analytic(rho::Real)::Float64
    # For a bivariate normal with correlation rho: I(X;Y) = -0.5 * log2(1 - rho^2)
    r = Float64(rho)
    return -0.5 * log2(1.0 - r^2)
end

function gaussian_mi_histogram_demo(; rho::Real=0.5, n_samples::Int=200_000, n_bins::Int=64, seed::Int=1234)
    rng = MersenneTwister(seed)

    # Generate correlated standard normal samples.
    x = randn(rng, n_samples)
    z = randn(rng, n_samples)
    y = rho .* x .+ sqrt(1 - rho^2) .* z

    mi_est = mutual_information_histogram(x, y; n_bins=n_bins)
    mi_true = gaussian_mi_analytic(rho)

    return (mi_est=mi_est, mi_true=mi_true)
end

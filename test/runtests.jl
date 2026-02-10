using Test

# Load the core modules directly (project-style repository, no Julia package required).
include(joinpath(@__DIR__, "..", "src", "core", "entanglement", "entanglement.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "random_unitaries", "random_unitaries.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "circuits", "circuits.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "utilities", "utilities.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "mutual_information", "mutual_information.jl"))

using .entanglement
using .random_unitaries
using .circuits
using .utilities
using .mutual_information

@testset "entanglement" begin
    include("test_entanglement_negativity.jl")
end

@testset "circuits" begin
    include("test_circuits_unitarity.jl")
end

@testset "utilities" begin
    include("test_utilities_numerical_integration.jl")
end

@testset "mutual_information" begin
    include("test_mutual_information.jl")
end

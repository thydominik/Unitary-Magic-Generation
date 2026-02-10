using Test

# Load the core modules directly (project-style repository, no Julia package required).
include(joinpath(@__DIR__, "..", "src", "core", "entanglement", "entanglement.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "random_unitaries", "random_unitaries.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "circuits", "circuits.jl"))

using .entanglement
using .random_unitaries
using .circuits

@testset "entanglement" begin
    include("test_entanglement_negativity.jl")
end

@testset "circuits" begin
    include("test_circuits_unitarity.jl")
end

using Test

# Load the core modules directly (project-style repository, no Julia package required).
include(joinpath(@__DIR__, "..", "src", "core", "entanglement", "entanglement.jl"))

using .entanglement

@testset "entanglement" begin
    include("test_entanglement_negativity.jl")
end

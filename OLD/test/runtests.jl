using Test
using UnitaryMagicGeneration

@testset "UnitaryMagicGeneration.jl" begin
    @testset "Core Module" begin
        include("core/test_magic.jl")
        include("core/test_random_unitaries.jl")
        include("core/test_entanglement.jl")
        include("core/test_entropy.jl")
    end
end

using Test
using LinearAlgebra
using Random

@testset "brickwall_unitary_matrix is unitary" begin
    rng = MersenneTwister(1234)

    # Use random_unitaries for deterministic 2-qubit gates.
    sampler = r -> random_unitaries.cue_matrix(4; rng=r)

    u2 = circuits.brickwall_unitary_matrix(2; depth=3, two_qubit_gate_sampler=sampler, rng=rng)
    @test size(u2) == (4, 4)
    @test isapprox(u2' * u2, Matrix{ComplexF64}(I, 4, 4); atol=1e-10, rtol=0)

    rng2 = MersenneTwister(5678)
    u3 = circuits.brickwall_unitary_matrix(3; depth=4, two_qubit_gate_sampler=sampler, rng=rng2)
    @test size(u3) == (8, 8)
    @test isapprox(u3' * u3, Matrix{ComplexF64}(I, 8, 8); atol=1e-10, rtol=0)
end

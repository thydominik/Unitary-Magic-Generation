using Test
using UnitaryMagicGeneration
using LinearAlgebra

@testset "RandomUnitaries Module" begin
    
    @testset "cue_matrix" begin
        # Test dimensions
        for dim in [2, 4, 8]
            U = cue_matrix(dim)
            @test size(U) == (dim, dim)
        end
        
        # Test unitarity
        U = cue_matrix(4)
        @test U * U' ≈ I atol=1e-10
        @test U' * U ≈ I atol=1e-10
        
        # Test determinant (should have |det(U)| = 1)
        @test abs(det(U)) ≈ 1.0 atol=1e-10
        
        # Test beta parameter (GOE, GUE, GSE)
        U_goe = cue_matrix(4, beta=1)
        U_gue = cue_matrix(4, beta=2)
        U_gse = cue_matrix(4, beta=4)
        
        @test U_goe * U_goe' ≈ I atol=1e-10
        @test U_gue * U_gue' ≈ I atol=1e-10
        @test U_gse * U_gse' ≈ I atol=1e-10
        
        # Test randomness (two calls should give different matrices)
        U1 = cue_matrix(4)
        U2 = cue_matrix(4)
        @test !(U1 ≈ U2)
        
        # Test invalid inputs
        @test_throws ArgumentError cue_matrix(0)
        @test_throws ArgumentError cue_matrix(-1)
    end
    
    @testset "generate_regular_unitary_circuit" begin
        # Test dimensions
        for n in [2, 3, 4]
            U = generate_regular_unitary_circuit(n)
            expected_dim = 2^n
            @test size(U) == (expected_dim, expected_dim)
        end
        
        # Test unitarity
        U = generate_regular_unitary_circuit(3)
        @test U * U' ≈ I atol=1e-10
        @test U' * U ≈ I atol=1e-10
        
        # Test determinant
        @test abs(det(U)) ≈ 1.0 atol=1e-10
        
        # Test randomness
        U1 = generate_regular_unitary_circuit(3)
        U2 = generate_regular_unitary_circuit(3)
        @test !(U1 ≈ U2)
        
        # Test invalid inputs
        @test_throws ArgumentError generate_regular_unitary_circuit(0)
        @test_throws ArgumentError generate_regular_unitary_circuit(-1)
    end
    
    @testset "generate_brickwall_unitary_circuit" begin
        # Test dimensions
        for n in [2, 3, 4]
            U = generate_brickwall_unitary_circuit(n)
            expected_dim = 2^n
            @test size(U) == (expected_dim, expected_dim)
        end
        
        # Test unitarity
        U = generate_brickwall_unitary_circuit(4)
        @test U * U' ≈ I atol=1e-10
        @test U' * U ≈ I atol=1e-10
        
        # Test determinant
        @test abs(det(U)) ≈ 1.0 atol=1e-10
        
        # Test custom depth
        U_short = generate_brickwall_unitary_circuit(4, depth=2)
        U_long = generate_brickwall_unitary_circuit(4, depth=20)
        
        @test U_short * U_short' ≈ I atol=1e-10
        @test U_long * U_long' ≈ I atol=1e-10
        
        # Test randomness
        U1 = generate_brickwall_unitary_circuit(3)
        U2 = generate_brickwall_unitary_circuit(3)
        @test !(U1 ≈ U2)
        
        # Test invalid inputs
        @test_throws ArgumentError generate_brickwall_unitary_circuit(0)
        @test_throws ArgumentError generate_brickwall_unitary_circuit(-1)
        @test_throws ArgumentError generate_brickwall_unitary_circuit(3, depth=0)
        @test_throws ArgumentError generate_brickwall_unitary_circuit(3, depth=-1)
    end
    
    @testset "Circuit properties" begin
        # Test that brickwall circuits preserve normalization
        n = 3
        U = generate_brickwall_unitary_circuit(n)
        
        # Apply to random state
        Random.seed!(42)
        ψ = randn(ComplexF64, 2^n)
        ψ /= norm(ψ)
        
        ψ_evolved = U * ψ
        
        # Check normalization preserved
        @test norm(ψ_evolved) ≈ 1.0 atol=1e-10
    end
end

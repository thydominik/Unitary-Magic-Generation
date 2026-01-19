using Test
using UnitaryMagicGeneration
using LinearAlgebra

@testset "Entanglement Module" begin
    
    @testset "reduced_density_matrix" begin
        # Test with 2-qubit product state |00⟩
        ψ = [1.0 + 0.0im, 0.0, 0.0, 0.0]
        
        # Trace out qubit 2
        ρ_reduced, kept = reduced_density_matrix(ψ, [2])
        @test kept == [1]
        @test size(ρ_reduced) == (2, 2)
        @test ρ_reduced ≈ [1.0 0.0; 0.0 0.0] atol=1e-10
        
        # Trace out qubit 1
        ρ_reduced, kept = reduced_density_matrix(ψ, [1])
        @test kept == [2]
        @test ρ_reduced ≈ [1.0 0.0; 0.0 0.0] atol=1e-10
        
        # Test with Bell state |Φ+⟩ = (|00⟩ + |11⟩)/√2
        ψ_bell = [1.0, 0.0, 0.0, 1.0] / sqrt(2)
        ρ_bell, _ = reduced_density_matrix(ψ_bell, [1])
        
        # Should be maximally mixed
        @test ρ_bell ≈ [0.5 0.0; 0.0 0.5] atol=1e-10
        
        # Test 3-qubit GHZ state |000⟩ + |111⟩
        ψ_ghz = zeros(ComplexF64, 8)
        ψ_ghz[1] = ψ_ghz[8] = 1 / sqrt(2)
        
        # Trace out qubit 3
        ρ_12, kept = reduced_density_matrix(ψ_ghz, [3])
        @test kept == [1, 2]
        @test size(ρ_12) == (4, 4)
        
        # Test invalid inputs
        @test_throws ArgumentError reduced_density_matrix([1.0, 0.0, 0.0], [1])  # Non-power-of-2
        @test_throws ArgumentError reduced_density_matrix(ψ, [3])  # Invalid qubit index
        @test_throws ArgumentError reduced_density_matrix(ψ, [1, 1])  # Duplicate indices
    end
    
    @testset "von_neumann_entropy" begin
        # Test with pure state density matrix
        ψ = [1.0 + 0.0im, 0.0]
        ρ_pure = ψ * ψ'
        S_pure = von_neumann_entropy(ρ_pure)
        @test S_pure ≈ 0.0 atol=1e-10
        
        # Test with maximally mixed state (1 qubit)
        ρ_mixed = [0.5 0.0; 0.0 0.5]
        S_mixed = von_neumann_entropy(ρ_mixed)
        @test S_mixed ≈ 1.0 atol=1e-10  # log2(2) = 1
        
        # Test with maximally mixed state (2 qubits)
        ρ_mixed_2 = Matrix(I, 4, 4) / 4
        S_mixed_2 = von_neumann_entropy(ρ_mixed_2)
        @test S_mixed_2 ≈ 2.0 atol=1e-10  # log2(4) = 2
        
        # Test with Bell state reduced density matrix (should be maximally mixed)
        ψ_bell = [1.0, 0.0, 0.0, 1.0] / sqrt(2)
        ρ_bell, _ = reduced_density_matrix(ψ_bell, [1])
        S_bell = von_neumann_entropy(ρ_bell)
        @test S_bell ≈ 1.0 atol=1e-10
    end
    
    @testset "calculate_entanglement" begin
        # Test with product state (no entanglement)
        ψ_product = [1.0 + 0.0im, 0.0, 0.0, 0.0]
        entropies, subsystems = calculate_entanglement(ψ_product; subsystems="All")
        
        @test length(entropies) == 1  # Only one cut for 2 qubits
        @test entropies[1] ≈ 0.0 atol=1e-10
        @test subsystems == [[2]]
        
        # Test with Bell state (maximal entanglement)
        ψ_bell = [1.0, 0.0, 0.0, 1.0] / sqrt(2)
        entropies_bell, _ = calculate_entanglement(ψ_bell; subsystems="All")
        
        @test entropies_bell[1] ≈ 1.0 atol=1e-10
        
        # Test with 3-qubit state
        ψ_3 = randn(ComplexF64, 8)
        ψ_3 /= norm(ψ_3)
        
        entropies_all, subsystems_all = calculate_entanglement(ψ_3; subsystems="All")
        @test length(entropies_all) == 2  # Two cuts: 1:2 and 2:1
        @test all(S >= 0.0 for S in entropies_all)
        
        # Test middle cut
        entropies_mid, subsystems_mid = calculate_entanglement(ψ_3; subsystems="Middle")
        @test length(entropies_mid) == 1
        @test entropies_mid[1] >= 0.0
        
        # Test invalid inputs
        @test_throws ArgumentError calculate_entanglement([1.0], subsystems="All")  # Only 1 qubit
        @test_throws ArgumentError calculate_entanglement(ψ_bell, subsystems="Invalid")
    end
    
    @testset "Entanglement properties" begin
        # Test that entanglement entropy is symmetric
        ψ = randn(ComplexF64, 8)
        ψ /= norm(ψ)
        
        # Trace out qubits 1,2 vs qubit 3
        ρ_3, _ = reduced_density_matrix(ψ, [1, 2])
        ρ_12, _ = reduced_density_matrix(ψ, [3])
        
        S_3 = von_neumann_entropy(ρ_3)
        S_12 = von_neumann_entropy(ρ_12)
        
        # For pure states, S(A) = S(B) for bipartition A:B
        @test S_3 ≈ S_12 atol=1e-10
    end
end

using Test
using UnitaryMagicGeneration
using LinearAlgebra

@testset "Entropy Module" begin
    
    @testset "renyi_entropy" begin
        # Test with pure state density matrix (α = 0)
        ψ = [1.0 + 0.0im, 0.0]
        ρ_pure = ψ * ψ'
        S0_pure = renyi_entropy(ρ_pure, 0.0)
        @test S0_pure ≈ 0.0 atol=1e-10  # log2(1) = 0 (rank = 1)
        
        # Test with maximally mixed state (1 qubit, various α)
        ρ_mixed = [0.5 0.0; 0.0 0.5]
        
        S0_mixed = renyi_entropy(ρ_mixed, 0.0)
        @test S0_mixed ≈ 1.0 atol=1e-10  # log2(2) = 1 (rank = 2)
        
        S1_mixed = renyi_entropy(ρ_mixed, 1.0)
        @test S1_mixed ≈ 1.0 atol=1e-10  # von Neumann limit
        
        S2_mixed = renyi_entropy(ρ_mixed, 2.0)
        @test S2_mixed ≈ 1.0 atol=1e-10  # Collision entropy
        
        # Test α = 1 gives von Neumann entropy
        ψ_bell = [1.0, 0.0, 0.0, 1.0] / sqrt(2)
        ρ_bell, _ = reduced_density_matrix(ψ_bell, [1])
        
        S1_bell = renyi_entropy(ρ_bell, 1.0)
        S_vn_bell = von_neumann_entropy(ρ_bell)
        @test S1_bell ≈ S_vn_bell atol=1e-10
        
        # Test that Rényi entropy is monotonic in α (for some cases)
        ρ_test = [0.7 0.0; 0.0 0.3]
        S_half = renyi_entropy(ρ_test, 0.5)
        S_one = renyi_entropy(ρ_test, 1.0)
        S_two = renyi_entropy(ρ_test, 2.0)
        
        # S_α decreases with increasing α for non-uniform distributions
        @test S_half >= S_one
        @test S_one >= S_two
        
        # Test invalid inputs
        @test_throws ArgumentError renyi_entropy(ρ_mixed, -1.0)  # Negative α
    end
    
    @testset "calculate_renyi_entropy" begin
        # Test with product state
        ψ_product = [1.0 + 0.0im, 0.0, 0.0, 0.0]
        S2_product = calculate_renyi_entropy(ψ_product, 2.0; subsystems="All")
        
        @test length(S2_product) == 1
        @test S2_product[1] ≈ 0.0 atol=1e-10
        
        # Test with Bell state
        ψ_bell = [1.0, 0.0, 0.0, 1.0] / sqrt(2)
        S2_bell = calculate_renyi_entropy(ψ_bell, 2.0; subsystems="All")
        
        @test S2_bell[1] >= 0.0
        
        # Compare α = 1 with von Neumann
        S1_bell = calculate_renyi_entropy(ψ_bell, 1.0; subsystems="All")
        S_vn_bell, _ = calculate_entanglement(ψ_bell; subsystems="All")
        
        @test S1_bell[1] ≈ S_vn_bell[1] atol=1e-10
        
        # Test with 3-qubit state (all cuts)
        ψ_3 = randn(ComplexF64, 8)
        ψ_3 /= norm(ψ_3)
        
        S2_all = calculate_renyi_entropy(ψ_3, 2.0; subsystems="All")
        @test length(S2_all) == 2
        @test all(S >= 0.0 for S in S2_all)
        
        # Test middle cut
        S2_mid = calculate_renyi_entropy(ψ_3, 2.0; subsystems="Middle")
        @test length(S2_mid) == 1
        @test S2_mid[1] >= 0.0
        
        # Test different α values
        for α in [0.0, 0.5, 1.0, 2.0, 5.0]
            Sα = calculate_renyi_entropy(ψ_bell, α; subsystems="All")
            @test all(S >= 0.0 for S in Sα)
        end
        
        # Test invalid inputs
        @test_throws ArgumentError calculate_renyi_entropy([1.0], 2.0, subsystems="All")
        @test_throws ArgumentError calculate_renyi_entropy(ψ_bell, -1.0, subsystems="All")
        @test_throws ArgumentError calculate_renyi_entropy(ψ_bell, 2.0, subsystems="Invalid")
    end
    
    @testset "Entropy relations" begin
        # Test that for pure bipartite states, Rényi entropies are equal across cuts
        ψ = randn(ComplexF64, 8)
        ψ /= norm(ψ)
        
        for α in [0.5, 1.0, 2.0]
            # Trace out 1 vs 2+3
            ρ_1, _ = reduced_density_matrix(ψ, [2, 3])
            ρ_23, _ = reduced_density_matrix(ψ, [1])
            
            Sα_1 = renyi_entropy(ρ_1, α)
            Sα_23 = renyi_entropy(ρ_23, α)
            
            # For pure states, entropies should match
            @test Sα_1 ≈ Sα_23 atol=1e-10
        end
    end
end

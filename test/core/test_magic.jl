using Test
using UnitaryMagicGeneration
using LinearAlgebra
using SparseArrays

@testset "Magic Module" begin
    
    @testset "pauli_matrix" begin
        # Test individual Pauli matrices
        @test pauli_matrix('I') ≈ [1 0; 0 1]
        @test pauli_matrix('X') ≈ [0 1; 1 0]
        @test pauli_matrix('Y') ≈ [0 -im; im 0]
        @test pauli_matrix('Z') ≈ [1 0; 0 -1]
        
        # Test invalid input
        @test_throws ArgumentError pauli_matrix('A')
    end
    
    @testset "generate_all_pauli_strings" begin
        # 1 qubit: I, X, Y, Z (4 strings)
        strings_1 = generate_all_pauli_strings(1)
        @test length(strings_1) == 4
        @test "I" in strings_1
        @test "X" in strings_1
        @test "Y" in strings_1
        @test "Z" in strings_1
        
        # 2 qubits: 16 strings
        strings_2 = generate_all_pauli_strings(2)
        @test length(strings_2) == 16
        @test "II" in strings_2
        @test "XX" in strings_2
        @test "YZ" in strings_2
        
        # Test invalid input
        @test_throws ArgumentError generate_all_pauli_strings(0)
        @test_throws ArgumentError generate_all_pauli_strings(-1)
    end
    
    @testset "pauli_operator_list" begin
        n = 2
        pauli_strings = generate_all_pauli_strings(n)
        pauli_ops = pauli_operator_list(pauli_strings, n)
        
        # Check correct number of operators
        @test length(pauli_ops) == 4^n
        
        # Check type
        @test all(op isa SparseMatrixCSC{ComplexF64} for op in pauli_ops)
        
        # Check dimensions
        @test all(size(op) == (2^n, 2^n) for op in pauli_ops)
        
        # Check specific operators
        # II should be identity
        ii_idx = findfirst(s -> s == "II", pauli_strings)
        @test Matrix(pauli_ops[ii_idx]) ≈ Matrix(I, 4, 4)
        
        # XX should be X ⊗ X
        xx_idx = findfirst(s -> s == "XX", pauli_strings)
        X = pauli_matrix('X')
        @test Matrix(pauli_ops[xx_idx]) ≈ kron(X, X)
        
        # Test invalid inputs
        @test_throws ArgumentError pauli_operator_list(["XY"], 3)  # Wrong n_qubits
        @test_throws ArgumentError pauli_operator_list(["XA"], 2)  # Invalid Pauli char
    end
    
    @testset "measure_magic_pure" begin
        n = 2
        pauli_strings = generate_all_pauli_strings(n)
        pauli_ops = pauli_operator_list(pauli_strings, n)
        
        # Test with computational basis state |00⟩
        ψ_00 = [1.0 + 0.0im, 0.0, 0.0, 0.0]
        magic_00, Ξ_00 = measure_magic_pure(ψ_00, pauli_ops)
        
        # Stabilizer state should have zero magic
        @test magic_00 ≈ 0.0 atol=1e-10
        @test length(Ξ_00) == length(pauli_ops)
        
        # Test with Bell state |Φ+⟩ = (|00⟩ + |11⟩)/√2
        ψ_bell = [1.0, 0.0, 0.0, 1.0] / sqrt(2)
        magic_bell, Ξ_bell = measure_magic_pure(ψ_bell, pauli_ops)
        
        # Bell state is also a stabilizer state
        @test magic_bell ≈ 0.0 atol=1e-10
        
        # Test with random state (should have non-zero magic in general)
        Random.seed!(42)
        ψ_random = randn(ComplexF64, 2^n)
        ψ_random /= norm(ψ_random)
        magic_random, Ξ_random = measure_magic_pure(ψ_random, pauli_ops)
        
        # Magic should be non-negative
        @test magic_random >= 0.0
        
        # Test invalid inputs
        @test_throws ArgumentError measure_magic_pure([1.0, 0.0], pauli_ops)  # Wrong dimension
        @test_throws ArgumentError measure_magic_pure(ψ_00, [])  # Empty operator list
    end
    
    @testset "measure_magic_mixed" begin
        n = 2
        pauli_strings = generate_all_pauli_strings(n)
        pauli_ops = pauli_operator_list(pauli_strings, n)
        
        # Test with pure state density matrix
        ψ = [1.0 + 0.0im, 0.0, 0.0, 0.0]
        ρ_pure = ψ * ψ'
        magic_pure, Ξ_pure = measure_magic_mixed(ρ_pure, pauli_ops)
        
        # Should match pure state result
        magic_pure_direct, _ = measure_magic_pure(ψ, pauli_ops)
        @test magic_pure ≈ magic_pure_direct atol=1e-10
        
        # Test with maximally mixed state
        ρ_mixed = Matrix(I, 2^n, 2^n) / 2^n
        magic_mixed, Ξ_mixed = measure_magic_mixed(ρ_mixed, pauli_ops)
        
        # Maximally mixed state should have specific magic value
        @test magic_mixed >= 0.0
        @test length(Ξ_mixed) == length(pauli_ops)
        
        # Test invalid inputs
        @test_throws ArgumentError measure_magic_mixed(zeros(2, 2), pauli_ops)  # Wrong dimension
    end
    
    @testset "Magic properties" begin
        n = 2
        pauli_strings = generate_all_pauli_strings(n)
        pauli_ops = pauli_operator_list(pauli_strings, n)
        
        # Test that Clifford states have zero magic
        clifford_states = [
            [1.0, 0.0, 0.0, 0.0],  # |00⟩
            [0.0, 1.0, 0.0, 0.0],  # |01⟩
            [0.0, 0.0, 1.0, 0.0],  # |10⟩
            [0.0, 0.0, 0.0, 1.0],  # |11⟩
            [1.0, 1.0, 0.0, 0.0] / sqrt(2),  # |+0⟩
            [1.0, 0.0, 0.0, 1.0] / sqrt(2),  # |Φ+⟩ Bell state
        ]
        
        for ψ in clifford_states
            magic, _ = measure_magic_pure(ComplexF64.(ψ), pauli_ops)
            @test magic ≈ 0.0 atol=1e-10
        end
    end
end

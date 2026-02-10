using Test
using LinearAlgebra

@testset "negativity (pure states)" begin
    # Bell state |Phi+> = (|00> + |11>)/sqrt(2)
    psi_bell = ComplexF64[1, 0, 0, 1] ./ sqrt(2)

    n = entanglement_negativity(psi_bell, [1])
    ln = logarithmic_negativity(psi_bell, [1])

    @test isapprox(n, 0.5; atol=1e-12, rtol=0)
    @test isapprox(ln, 1.0; atol=1e-12, rtol=0)

    # Product state |00> has zero entanglement
    psi_prod = ComplexF64[1, 0, 0, 0]

    n0 = entanglement_negativity(psi_prod, [1])
    ln0 = logarithmic_negativity(psi_prod, [1])

    @test isapprox(n0, 0.0; atol=1e-12, rtol=0)
    @test isapprox(ln0, 0.0; atol=1e-12, rtol=0)
end

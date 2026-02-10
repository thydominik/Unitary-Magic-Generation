"""
    random_unitaries

Core routines for generating random unitary matrices.

This module intentionally contains only unitary matrix sampling utilities.
Circuit construction lives in `src/core/circuits`.

All identifiers and docstrings are ASCII-only.

Implemented distributions
- Haar-random unitary matrices (CUE / U(N)) via the standard Ginibre-QR method.

Notes
- The Haar sampler is deterministic with respect to the provided RNG.
- No external dependencies are required.
"""
module random_unitaries

using LinearAlgebra
using Random

export cue_matrix,
    random_unitary_matrix

"""
    cue_matrix(dim; rng=Random.default_rng()) -> Matrix{ComplexF64}

Draw a Haar-random unitary matrix from U(dim).

Parameters
- `dim`: Matrix dimension (must be positive).
- `rng`: Random-number generator.

Returns
- `U`: A `dim x dim` unitary matrix.

Implementation
Uses the Ginibre-QR construction for Haar sampling on U(N):
1) Sample a complex Gaussian matrix Z.
2) Compute QR factorization Z = Q R.
3) Fix phases using diag(R) so the distribution is Haar.
"""
function cue_matrix(dim::Integer; rng::AbstractRNG=Random.default_rng())::Matrix{ComplexF64}
    dim > 0 || throw(ArgumentError("dim must be positive, got $dim"))

    z_re = randn(rng, dim, dim)
    z_im = randn(rng, dim, dim)
    z = ComplexF64.(z_re, z_im)

    f = qr(z)
    q = Matrix(f.Q)
    r = f.R

    d = diag(r)
    ph = similar(d)
    for i in eachindex(d)
        a = abs(d[i])
        ph[i] = a == 0 ? (1.0 + 0im) : (d[i] / a)
    end

    return q * Diagonal(ph)
end

"""
    random_unitary_matrix(n_qubits; rng=Random.default_rng()) -> Matrix{ComplexF64}

Draw a Haar-random unitary on `n_qubits` qubits, i.e., a matrix of size `2^n_qubits x 2^n_qubits`.
"""
function random_unitary_matrix(n_qubits::Integer; rng::AbstractRNG=Random.default_rng())::Matrix{ComplexF64}
    n_qubits > 0 || throw(ArgumentError("n_qubits must be positive, got $n_qubits"))
    return cue_matrix(2^n_qubits; rng=rng)
end

end # module random_unitaries

using Test
using Random

@testset "mutual_information_histogram basic" begin
    rng = MersenneTwister(123)

    n = 50_000
    x = randn(rng, n)

    # Dependent case.
    y_dep = x

    # Independent case (shuffle).
    y_ind = copy(x)
    shuffle!(rng, y_ind)

    mi_dep = mutual_information_histogram(x, y_dep; n_bins=40)
    mi_ind = mutual_information_histogram(x, y_ind; n_bins=40)

    @test mi_dep >= 0.0
    @test mi_ind >= 0.0

    # Dependent should be larger by a healthy margin.
    @test mi_dep > mi_ind + 0.2
end

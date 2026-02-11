using Test

@testset "double_integral_trapz" begin
    # Integrate f(x,y) = x + y over [0,1]x[0,1].
    # Exact result: integral_0^1 integral_0^1 (x+y) dx dy = 1.0.
    x = collect(range(0.0, 1.0; length=101))
    y = collect(range(0.0, 1.0; length=121))

    f = [x[i] + y[j] for i in eachindex(x), j in eachindex(y)]

    val = utilities.double_integral_trapz(f, x, y)
    @test isapprox(val, 1.0; atol=1e-6, rtol=0)
end

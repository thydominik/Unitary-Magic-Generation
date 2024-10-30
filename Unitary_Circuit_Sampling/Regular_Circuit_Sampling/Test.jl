t1 = time()
t2 = time()
t2 - t1
Pauli_Matrices = ["I", "X", "Y", "Z"]
n = 2

function timing(n)
    t1 = time()
    a = [join(s) for s in Iterators.product(fill(Pauli_Matrices, n)...)];
    t2 = time()

    return a, (t2 - t1)
end

a, dt = timing(12);
println(dt)

##

include("Magic.jl")
using .Measure_Magic
using LinearAlgebra
using Plots
N = 2
S = Measure_Magic.GenerateAllPauliStrings(N);

Magic = Vector()
t = Vector()
push!(t, time())
for i in 1:10000
    println(i)
    p = rand(2^N) + im * rand(2^N);
    p = p ./ norm(p);

    push!(Magic, Measure_Magic.MeasureMagic(p, S, 2));
    push!(t, time());
end

plot(sort(Magic)/log(5/2))
plot(sort(diff(t))[1:8])
a = Measure_Magic.PauliMatrix(3)

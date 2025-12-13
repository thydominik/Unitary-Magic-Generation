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
N = 8
S = Measure_Magic.GenerateAllPauliStrings(N);

Magic = Vector()
t = Vector()
push!(t, time())
for i in 1:1
    println(i)
    p = rand(2^N) + im * rand(2^N);
    p = p ./ norm(p);

    push!(Magic, Measure_Magic.MeasureMagic(p, S, 2));
    push!(t, time());
end

plot(sort(Magic)/log(2^N /2))
plot(sort(diff(t))[1:10])
a = Measure_Magic.PauliMatrix(3)



include("Magic.jl")
include("Random_Unitaries.jl")
using .Measure_Magic
using .Random_Unitary_Generation
using LinearAlgebra
using BenchmarkTools

@btime Random_Unitary_Generation.Generate_BW_Unitary_Ciruit(10, 1);


@btime Measure_Magic.GenerateAllPauliStrings(7);
print("\033c")
N = 8
S = Measure_Magic.GenerateAllPauliStrings(N);

PO = Measure_Magic.PauliOperatorList(S, N)
Base.summarysize(PO)/(1024^2)
Measure_Magic.MeasureMagic(rand(2^N), PO, 2)

A = rand(ComplexF64, 100, 100)  # 100x100 matrix of complex numbers
println("Data size (bytes): ", sizeof(A))
println("Total memory size (bytes): ", Base.summarysize(A))
Bas


using Random
using ProgressBars

include("Random_Unitaries.jl")
include("Magic.jl")
using .Random_Unitary_Generation
using .Measure_Magic

# Setting the seed for the random number generation
Seed = 1
Random.seed!(Seed)

# Sampling parameters
No_Samples = 2^20

No_Qubits = 1

Psi_0 = 1/sqrt(2^No_Qubits) * ones(2^No_Qubits);


Strings = Measure_Magic.GenerateAllPauliStrings(No_Qubits)
PauliOperators = Measure_Magic.PauliOperatorList(Strings, No_Qubits)
State = Psi_0
No_Qubits = Int(log2(length(State)))
HilbertSpaceDimension = 2^No_Qubits
α = 2
Ξ = Vector{Float64}(undef, 4^No_Qubits)
for PauliIndex in 1:4^No_Qubits
    Operator = PauliOperators[PauliIndex]
    println(Operator)
    println(Operator * State[:])
    println(conj(transpose(State[:])) * Operator * State[:])
    println((conj(transpose(State[:])) * Operator * State[:])^2)
    println(real((1/HilbertSpaceDimension) * (conj(transpose(State[:])) * Operator * State[:]))^2)
    Ξ[PauliIndex] = real((1/HilbertSpaceDimension) * (conj(transpose(State[:])) * Operator * State[:])^2)
end # FOR PauliIndex
Magic = (1 - α)^(-1) * log(sum(Ξ.^α)) - log(HilbertSpaceDimension)




using JLD2

a = 2
b = 3
c = 4
vector = rand(2^10)

fname = "RegularUnitaryCircuitMagic_$(a)_iter_$(b).jdl2"
@save fname vector a b c

data = JLD2.load(fname) 


data = JLD2.load("RegularUnitaryCircuitMagicSampled_N_6_Samples_1048576_Seed_1.jdl2")
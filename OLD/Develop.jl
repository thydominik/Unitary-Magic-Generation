###
# Develop is an actively changed, non-essential part of the code.
# It is for me to play around with functions, tests, and check datafiles before usage.
# Anything found here is not part of the core library and cannot be taken seriously, not to mention these parts must be prone to bugs.
###

using LinearAlgebra
using Trace
const P = [
    [1 0; 0 1],  # Identity
    [0 1; 1 0],  # Pauli X
    [0 -im; im 0], # Pauli Y
    [1 0; 0 -1]  # Pauli Z
]

Psi1 = rand(ComplexF64, 2); Psi1 /= norm(Psi1);
Psi2 = rand(ComplexF64, 2); Psi2 /= norm(Psi2);

Psi = kron(Psi1, Psi2);

m1 = 1/4 * (conj(transpose(Psi[:])) * kron(P[1], I(2)) * Psi[:])^2
m2 = 1/4 * (conj(transpose(Psi[:])) * kron(P[2], I(2)) * Psi[:])^2
m3 = 1/4 * (conj(transpose(Psi[:])) * kron(P[3], I(2)) * Psi[:])^2
m4 = 1/4 * (conj(transpose(Psi[:])) * kron(P[4], I(2)) * Psi[:])^2

rho = Psi[:] * conj(transpose(Psi[:]))
j1 = 1/4 * sum(diag(kron(P[1], I(2)) * rho))^2
j2 = 1/4 * sum(diag(kron(P[2], I(2)) * rho))^2
j3 = 1/4 * sum(diag(kron(P[3], I(2)) * rho))^2
j4 = 1/4 * sum(diag(kron(P[4], I(2)) * rho))^2

M1 = -log(m1^2 + m2^2 + m3^2 + m4^2) - log(4)
M2 = -log(j1^2 + j2^2 + j3^2 + j4^2) - log(4)

current_dir = @__DIR__

# Unitary matrices
filepath = joinpath("Modules", "Magic.jl")
include(filepath)
# Magic
filepath = joinpath("Modules", "Random_Unitaries.jl")
include(filepath)
# Measure_Entanglement
filepath = joinpath("Modules", "Entanglement.jl")
include(filepath)

using .Measure_Entanglement
using .Measure_Magic
using .Random_Unitary_Generation
Strings = Dict()
Strings[1] = Measure_Magic.GenerateAllPauliStrings(1)
Strings[2] = Measure_Magic.GenerateAllPauliStrings(2)
PauliOperators = Dict()
PauliOperators[1] = Measure_Magic.PauliOperatorList(Strings[1], 1)
PauliOperators[2] = Measure_Magic.PauliOperatorList(Strings[2], 2)

Measure_Magic.MeasureMagic(Psi1, PauliOperators, 2)


rho1 = Psi1 * Psi1'
rho2 = Psi2 * Psi2'
rho = 

Measure_Entanglement.reduced_density_matrix(kron(Psi1, Psi2), collect(combinations(1:1, 1)))
using Combinatorics

reduced_matrix = Measure_Entanglement.reduced_density_matrix(Psi, [1])

(reduced_matrix[1] - rho1)

Measure_Magic.MeasureMagic_Mixed(reduced_matrix[1], PauliOperators[1])
Measure_Magic.MeasureMagic_Pure(Psi1, PauliOperators[1])


using LinearAlgebra
N = 2

Strings             = Dict()
Strings[1]          = Measure_Magic.GenerateAllPauliStrings(N)
PauliOperators[1]   = Measure_Magic.PauliOperatorList(Strings[1], N)
d = 2^N
q = Dict()
for i in 1:4^N
    q[i] = kron(fill(PauliOperators[1][i], 4)...)
end
Q = zeros(d^4, d^4)
for i in 1:4^N
    Q += 1/d^2 * q[i]
end



Magic = Dict()

psi = [1; 1; 1; 1];
psi /= norm(psi)
rho = psi * psi'
rho[1, 2]
rho[2, 1]
M = 100
E1 = LinRange(rho[1, 2], 0, M)
E2 = LinRange(rho[1, 3], 0, M)
E3 = LinRange(rho[1, 4], 0, M)
E4 = LinRange(rho[2, 3], 0, M)
E5 = LinRange(rho[2, 4], 0, M)
E6 = LinRange(rho[3, 4], 0, M)

MQ = Vector{Float64}()
MM = Vector{Float64}()
MP = Vector{Float64}()
for i in 1:M
    rho[1, 2] = rho[2, 1] = E1[i]
    rho[1, 3] = rho[3, 1]= E2[i]
    rho[1, 4] = rho[4, 1]= E3[i]
    rho[2, 3] = rho[3, 2]= E4[i]
    rho[2, 4] = rho[4, 2]= E5[i]
    rho[3, 4] = rho[4, 3]= E6[i]
    println(rho)
    R = kron(fill(rho, 4)...)
    push!(MQ, -log2(d * sum(diag(Q * R))) - (-log2(sum(eigvals(rho).^2))))
    push!(MM, Measure_Magic.MeasureMagic_Mixed(rho, PauliOperators[1]))
    push!(MP, Measure_Magic.MeasureMagic_Pure(psi, PauliOperators[1]))
end


using Plots

plot(MQ ./ log2(3/2))
scatter!(MM ./ log2(3/2))
plot!(MP ./ log2(3/2))
plot!(ylims=[0, 1])


using JLD2

data = JLD2.load("D:\\Data\\Random_Unitary_Magic_Generation\\N5\\Regular\\RegularUnitaryCircuitMagicSampled_N_5_Samples_1048576_Seed_1.jld2")

histogram(data["Magic"])






D:\Data\Random_Unitary_Magic_Generation\N8

Data_1 = JLD2.load("D:\\Data\\Random_Unitary_Magic_Generation\\N8\\RegularUnitaryCircuitMagicSampled_N_8_Samples_425984.jld2")

unique(Data_1["Svn"])



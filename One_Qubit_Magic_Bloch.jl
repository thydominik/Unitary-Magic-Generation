using Random
using ProgressBars
using JLD2
using Base.Threads
using MAT
using ProgressBars

current_dir = @__DIR__

# Unitary matrices
filepath = joinpath(current_dir, "..", "Unitary-Magic-Generation/Modules", "Magic.jl")
include(filepath)

# Magic
filepath = joinpath(current_dir, "..", "Unitary-Magic-Generation/Modules", "Random_Unitaries.jl")
include(filepath)
# Measure_Entanglement
filepath = joinpath(current_dir, "..", "Unitary-Magic-Generation/Modules", "Entanglement.jl")
include(filepath)

using .Measure_Entanglement
using .Measure_Magic
using .Random_Unitary_Generation

Seed = 1
# Setting the seed for the random number generation
Random.seed!(Seed)
# Set the number of qubits first
No_Qubits = 1

Strings = Measure_Magic.GenerateAllPauliStrings(No_Qubits)
PauliOperators = Measure_Magic.PauliOperatorList(Strings, No_Qubits)

No_Samples = 2^13
theta = range(0, pi, Int(No_Samples / 2))
phi     = range(0, 2 * pi, Int(No_Samples) )
Magic   = zeros(length(theta), length(phi))

for i in ProgressBar(1:length(theta))
    for j in 1:length(phi)
        psi = [cos(theta[i] / 2), sin(theta[i] / 2) * exp(im * phi[j])]
        Magic[i, j] = Measure_Magic.MeasureMagic_Pure(psi, PauliOperators, 2)[1]
    end
end

using Plots
@save joinpath(current_dir, "Magic_Bloch_Data.jld2") Magic theta phi
matwrite(joinpath(current_dir, "Magic_Bloch_Data.mat"), Dict(
    "Magic" => Magic,
    "theta" => collect(theta),
    "phi" => collect(phi)
))


using Plots
plotlyjs()

# Create meshgrid for θ and ϕ
Θ = [θ for θ in theta, φ in phi]  # (length(theta), length(phi))
Φ = [φ for θ in theta, φ in phi]

# Convert spherical to Cartesian coordinates
X = sin.(Θ) .* cos.(Φ)
Y = sin.(Θ) .* sin.(Φ)
Z = cos.(Θ)

# Flatten everything for scatter plot
x = vec(X)
y = vec(Y)
z = vec(Z)
colors = vec(Magic)  # Flatten Magic values to match points

# Plot as a 3D scatter plot
scatter3d(x, y, z, marker_z=colors, c=:viridis, ms=4, legend=false,
          xlabel="X", ylabel="Y", zlabel="Z", title="Magic Scatter on Bloch Sphere")


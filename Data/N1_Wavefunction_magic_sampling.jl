using Random
using ProgressBars
using JLD2
using Base.Threads

current_dir = @__DIR__

# Unitary matrices
filepath = joinpath(current_dir, "Modules", "Magic.jl")
include(filepath)
# Magic
filepath = joinpath(current_dir, "Modules", "Random_Unitaries.jl")
include(filepath)
# Measure_Entanglement
filepath = joinpath(current_dir, "Modules", "Entanglement.jl")
include(filepath)

using .Random_Unitary_Generation
using .Measure_Magic
using .Measure_Entanglement

# Ψ =  (cos(θ₁)|0⟩ + e^(iφ₁)sin(θ₁)|1⟩) ⊗ (cos(θ₂)|0⟩ + e^(i(φ₂+δ))sin(θ₂)|1⟩)

function one_qubit_wavefunction(theta1, phi1)
    # Define basis states
    zero = [1.0 + 0.0im, 0.0 + 0.0im]
    one = [0.0 + 0.0im, 1.0 + 0.0im]
    
    # Parameterize the state
    state1 = cos(theta1) * zero + exp(1im * phi1) * sin(theta1) * one
    
    return state1
end

m = 2000


# Create sample points for each parameter
theta1_samples  = range(0, π,   length=m)
phi1_samples    = range(0, 2π,  length=2*m)

Strings         = Measure_Magic.GenerateAllPauliStrings(1)
PauliOperators  = Measure_Magic.PauliOperatorList(Strings, 1)

Magic = Vector{Float64}()
for theta1 in ProgressBar(theta1_samples), phi1 in phi1_samples
    wf = one_qubit_wavefunction(theta1, phi1)
    push!(Magic, Measure_Magic.MeasureMagic(wf, PauliOperators))
end

histogram(Magic)
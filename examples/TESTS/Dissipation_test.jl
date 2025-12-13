using LinearAlgebra
using SparseArrays
using Random
using ProgressBars
using JLD2
using Base.Threads
using MAT
using StatsBase: sample
using Plots

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
begin
        
    
    t = LinRange(0, 10, 500)
    α = 2
    γ = 1;
    Ψ = rand(2); Ψ /= norm(Ψ)
    
    
    Strings = Measure_Magic.GenerateAllPauliStrings(1)
    PauliOperators = Measure_Magic.PauliOperatorList(Strings, 1)
    ρ = Ψ * Ψ'
    Magic = Vector{Float64}()
    for i in 1:length(t)
        dissipation = exp(-4 * t[i])
        temp_rho = Ψ * Ψ'
        temp_rho[1, 2] *= dissipation
        temp_rho[2, 1] *= dissipation
        push!(Magic, Measure_Magic.MeasureMagic_Mixed(temp_rho, PauliOperators, α))
    end
end
plot(t, Magic, label = "Magic", xlabel = "Time", ylabel = "Magic", title = "Dissipation Test", ylims = (0, 1.1))
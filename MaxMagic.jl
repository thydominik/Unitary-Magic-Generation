using LinearAlgebra
using SparseArrays
using Random
using ProgressBars
using JLD2
using Base.Threads
using MAT
using StatsBase: sample

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

#N = 3 maximum magic state
Mmax = 0
States = Vector()
Strings = Measure_Magic.GenerateAllPauliStrings(4);
PauliOperators = Measure_Magic.PauliOperatorList(Strings, 4);
a = 2;
for c1 in [1, im, 1 + im, 0]
    for c2 in [1, im, 1+im, 0]
        for c3 in [1, im, 1+im, 0]
            for c4 in [1, im, 1+im, 0]
                for c5 in [1, im, 1+im, 0]
                    for c6 in [1, im, 1+im, 0]
                        for c7 in [1, im, 1+im, 0]
                            for c8 in [1, im, 1+im, 0]
                                for c9 in [1, im, 1 + im, 0]
                                    for c10 in [1, im, 1+im, 0]
                                        for c11 in [1, im, 1+im, 0]
                                            for c12 in [1, im, 1+im, 0]
                                                for c13 in [1, im, 1+im, 0]
                                                    for c14 in [1, im, 1+im, 0]
                                                        for c15 in [1, im, 1+im, 0]
                                                            for c16 in [1, im, 1+im, 0]
                                                                State = [c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16]
                                                                State /= norm(State)
                                                                M = Measure_Magic.MeasureMagic_Pure(State, PauliOperators, a)
                                                                if M > log2(17/2) - 0.1
                                                                    Mmax = M
                                                                    push!(States, State)
                                                                    println(M)
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

                                



vals = Vector()
for i in 1:length(PauliOperators)
    push!(vals, real.(conj(States[end]) * PauliOperators[i] * Transpose(States[end]))^2)
end

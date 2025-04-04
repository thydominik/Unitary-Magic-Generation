# Sampling of N qubit random unitary Brick Wall circuits.
# There are two paramters to a BW circuit
# D - Depth of the circuit
# N - Number of qubits
# 
# The starting state is irrelevant, but due to the CUE matrices used, any state suffices, I've chosen a 0 magic state.
# |Ψ⟩ = ⊗_i^N |+>
# U|Ψ⟩ = |Haar State⟩
begin
    using Random
    using ProgressBars
    using JLD2
    using Base.Threads
end

begin
    current_dir = @__DIR__

    # Unitary matrices
    filepath = joinpath(current_dir, "..", "Modules", "Magic.jl")
    include(filepath)
    # Magic
    filepath = joinpath(current_dir, "..", "Modules", "Random_Unitaries.jl")
    include(filepath)
    #include("C:\\Dominik\\PhD\\Projects\\Unitary-Magic-Generation\\Modules\\Magic.jl")
    # Measure_Entanglement
    filepath = joinpath(current_dir, "..", "Modules", "Entanglement.jl")
    include(filepath)

    filepath = joinpath(current_dir, "..", "Modules", "Entropy.jl")
    include(filepath)
    using .Measure_Entanglement
    using .Measure_Magic
    using .Random_Unitary_Generation
    using .Measure_Entropy

end

function BrickWall_Circuit_Magic_sampling(No_Qubits::Int, Depth::Int, No_Samples::Int, Seed::Int)
    # Setting the seed for the random number generation
    Random.seed!(Seed)

    # Sampling parameters
    No_Samples = No_Samples

    # Set the number of qubits first
    No_Qubits = No_Qubits

    Depth = Depth
 
    Psi_0 = 1/sqrt(2^No_Qubits) * ones(2^No_Qubits);

    Strings = Measure_Magic.GenerateAllPauliStrings(No_Qubits)
    PauliOperators = Measure_Magic.PauliOperatorList(Strings, No_Qubits)

    Magic1          = Vector{Float64}()
    Magic           = Vector{Vector{Float64}}()
    #Xi              = Vector{Vector{Float64}}()
    Entanglement    = Vector{Vector{Float64}}()
    SubSystems      = Vector{Vector{Vector{Int}}}()
    REntropy        = Vector{Vector{Float64}}()

    for i in ProgressBar(1:No_Samples)    
        U       = Random_Unitary_Generation.Generate_BW_Unitary_Circuit(No_Qubits, Depth);
        State   = U * Psi_0
        
        # Von Neumann Netanglement entropy
        SVN, SubSys = Measure_Entanglement.calculate_entanglement(State, subsystems = "All")
        # Different alpha number renyi entropy calculations
        RE1 = Measure_Entropy.calculate_Renyi_entropy(State, 1.0001)
        RE1_5 = Measure_Entropy.calculate_Renyi_entropy(State, 1.5)
        RE2 = Measure_Entropy.calculate_Renyi_entropy(State, 2.0)
        RE3 = Measure_Entropy.calculate_Renyi_entropy(State, 3.0)
        RE4 = Measure_Entropy.calculate_Renyi_entropy(State, 4.0)
        # Magic calculations
        M1, ~ = Measure_Magic.MeasureMagic_Pure(State, PauliOperators, 1)
        M, ~  = Measure_Magic.MeasureMagic_Pure(State, PauliOperators, [1.5, 2, 3, 4])
        
        push!(Magic1, M1)
        push!(Magic, M)
        # push!(Xi, Ξ)
        push!(Entanglement, SVN)
        if i == 1
            push!(SubSystems, SubSys)
        end
        push!(REntropy, [RE1; RE1_5; RE2; RE3; RE4])

    end

    fname = "BWUnitaryCircuit_Magic_Entanglement_Entropy_Sampled_N_$(No_Qubits)_D_$(Depth)_Samples_$(No_Samples)_Seed_$(Seed).jld2"
    @save fname Magic1 Magic Entanglement SubSystems REntropy Psi_0 No_Samples No_Qubits Depth Seed
end


N = 5
for D in 1:25
    Partitions = 1
    div = Int(log2(Partitions))
    for s in 1:Partitions
        BrickWall_Circuit_Magic_sampling(N, D, 2^19, s)
    end
end

a = load("BWUnitaryCircuit_Magic_Entanglement_Entropy_Sampled_N_5_D_1_Samples_1024_Seed_1.jld2")


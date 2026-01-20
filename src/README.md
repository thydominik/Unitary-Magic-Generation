# Source Code Organization

This document describes the organization of the `src/` directory.

## Active Module Structure

The codebase is organized into functional modules:

### `src/core/` - Core Quantum Computing Functionality

Core algorithms for quantum state analysis:

- **`Core.jl`**: Main module that aggregates all core submodules
- **`Magic.jl`**: Magic state quantification using stabilizer Rényi entropy
  - `pauli_operator_list()`: Build tensor-product Pauli operators
  - `measure_magic_pure()`: Magic quantification for pure states
  - `measure_magic_mixed()`: Magic quantification for mixed states
- **`Entanglement.jl`**: Entanglement quantification
  - `reduced_density_matrix()`: Partial trace operations
  - `von_neumann_entropy()`: Entropy calculation
  - `calculate_entanglement()`: Multi-cut entanglement analysis
- **`Entropy.jl`**: Various entropy measures
  - `shannon_entropy()`: Shannon entropy in bits
  - `renyi_entropy()`: Rényi entropy (generalized)
- **`RandomUnitaries.jl`**: Random unitary generation
  - `random_unitary_haar()`: Haar-random unitaries
  - `random_state()`: Random quantum states

### `src/analysis/` - Statistical Analysis Tools

Information-theoretic and statistical analysis:

- **`Analysis.jl`**: Main analysis module
- **`mutual_information.jl`**: Mutual information computation
  - `MI()`: Fixed-range mutual information (assumes data in [0,10])
  - `MIn()`: Data-driven mutual information (normalized to data range)

### `src/utils/` - Utility Functions

General-purpose numerical tools:

- **`Utils.jl`**: Main utils module
- **`numerical_integration.jl`**: Numerical integration routines
  - `double_integral_trapz()`: 2D trapezoidal rule integration

### `src/circuits/` - Circuit Operations (Placeholder)

Circuit generation and manipulation (to be populated):

- **`Circuits.jl`**: Placeholder for future circuit functionality

## Main Entry Point

**`src/UnitaryMagicGeneration.jl`**: Top-level package module that:
- Loads all submodules in correct dependency order
- Re-exports public APIs
- Provides unified interface

## Legacy Code (Archived)

Old code has been moved to `../archive/` and should NOT be used:

- ❌ `src/Analysis/` (capitalized) - Old version with data loading at init
- ❌ `src/analytics/` - Duplicate MI implementations
- ❌ `src/Mutual_information/` - Old standalone scripts
- ❌ `src/Modules/` - Legacy module structure
- ❌ Standalone `.jl` files in `src/` root

## Usage Examples

### Using Core Functionality

```julia
using UnitaryMagicGeneration
using LinearAlgebra

# Generate random quantum state
n_qubits = 3
ψ = random_state(n_qubits)

# Compute magic
pauli_ops = pauli_operator_list(generate_all_pauli_strings(n_qubits), n_qubits)
magic, probs = measure_magic_pure(ψ, pauli_ops)

# Compute entanglement
entropies, subsystems = calculate_entanglement(ψ; subsystems="All")
```

### Using Analysis Tools

```julia
using UnitaryMagicGeneration

# Create correlated data
x = randn(10000) .+ 5
y = 0.7 .* x .+ 0.3 .* randn(10000)

# Compute mutual information
mi, = MI(x, y, 2^10)  # Fixed range [0,10]
println("Mutual Information: $mi bits")

mi_norm, = MIn(x, y, 2^10)  # Data-driven range
println("Normalized MI: $mi_norm bits")
```

## Development Guidelines

1. **All new code goes in the active modules** (`core/`, `analysis/`, `utils/`, `circuits/`)
2. **Never add code to `src/` root** - use appropriate module directory
3. **Follow naming conventions**:
   - Module files: `PascalCase.jl`
   - Function names: `snake_case()`
   - Type names: `PascalCase`
4. **Add comprehensive docstrings** to all exported functions
5. **Write tests** in `test/` directory for new functionality
6. **Update this README** when adding new modules or significant functionality

## Module Dependencies

Dependency order (earlier modules don't depend on later ones):

1. `utils/` - No internal dependencies
2. `core/` - Depends on `utils/`
3. `analysis/` - Depends on `utils/`
4. `circuits/` - May depend on `core/` when implemented

## Questions?

See the main repository README or individual module docstrings for more details.

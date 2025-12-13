# Refactoring Guide: Modular Code Structure

## Overview

This refactoring reorganizes the codebase into a modular, maintainable structure following Julia best practices:

- **Separation of Concerns**: Each functionality lives in its own module
- **Clear Dependencies**: Explicit module imports and dependencies
- **Improved Testability**: Isolated functions are easier to test
- **Better Documentation**: Comprehensive docstrings with examples
- **Type Safety**: Explicit type hints for better error messages

## Directory Structure

```
Unitary-Magic-Generation/
├─ src/                           # Source code
│  ├─ UnitaryMagic.jl              # Main module entry point
│  ├─ utils/                      # Utility functions
│  │  └─ numerical_integration.jl  # 2D trapezoidal integration
│  └─ analytics/                   # Analysis modules  
│     └─ mutual_information.jl     # MI computation & analysis
│
├─ examples/                      # Example scripts
│  └─ mutual_information_analysis.jl  # MI analysis pipeline
│
├─ test/                         # Test suite (future)
│  └─ test_mutual_information.jl  # MI function tests (future)
│
├─ Modules/                      # Existing modules (legacy)
├─ REFACTORING_GUIDE.md           # This file
└─ README.md                      # Project README
```

## Module Architecture

### 1. `NumericalIntegration` (`src/utils/numerical_integration.jl`)

**Purpose**: Provides numerical integration utilities for scientific computing.

**Exported Functions**:
- `double_integral_trapz(f, x, y)`: 2D trapezoidal rule integration

**Key Features**:
- Type-safe input validation
- Comprehensive error messages
- Well-documented with docstrings

**Usage**:
```julia
using UnitaryMagic
f = randn(100, 100)
x = range(0, 1, 100)
y = range(0, 1, 100)
result = double_integral_trapz(f, x, y)
```

### 2. `MutualInformationAnalysis` (`src/analytics/mutual_information.jl`)

**Purpose**: Compute mutual information between two distributions using histogram-based estimation.

**Exported Functions**:
- `MI(x, y, N_bins)`: Fixed-range MI (0 to 10)
- `MIn(x, y, N_bins)`: Data-driven MI (normalized)

**Return Values**:
Both functions return a tuple:
```julia n
(mi, px_individual, px_marginal, py_individual, py_marginal, pxy, Integrand)
```

**Theory**:
Mutual Information is computed as:
```
MI(X;Y) = ∫∫ p(x,y) log₂(p(x,y) / (p(x)p(y))) dx dy
```

**Usage**:
```julia
using UnitaryMagic
x = randn(1000) .+ 5  # Magic distribution
y = randn(1000) .+ 5  # Entanglement distribution

# Fixed range [0, 10]
mi, px_ind, px_marg, py_ind, py_marg, pxy, integrand = MI(x, y, 2^12)

# Data-driven range
mi_norm, = MIn(x, y, 2^12)
```

### 3. `UnitaryMagic` (`src/UnitaryMagic.jl`)

**Purpose**: Main package module that aggregates all functionality.

**Exported Functions**:
- All functions from submodules: `MI`, `MIn`, `double_integral_trapz`

**Usage**:
```julia
using UnitaryMagic

# All functions available directly
result = MI(x, y, 100)
```

## Key Improvements from Original Code

### 1. **Code Organization**

**Before**: All functions mixed in one large file
```julia
# MutualInformation.jl (7531 bytes)
function double_integral_trapz(...)
end
function MI(...)
end  
function MIn(...)
end
# + 100+ lines of example/analysis code
```

**After**: Separated by functionality
```
src/utils/numerical_integration.jl     (module with 1 function)
src/analytics/mutual_information.jl    (module with 2 functions)  
examples/mutual_information_analysis.jl (analysis pipeline)
```

### 2. **Documentation**

**Before**: Minimal or no comments
```julia
function MI(x, y, N_bins)
    m = length(x)
    # .. code without explanation ..
end
```

**After**: Comprehensive docstrings
```julia
"""
    MI(x::Vector, y::Vector, N_bins::Int)::Tuple

Compute mutual information using fixed range [0, 10]...

# Arguments
- `x::Vector`: First data array
- `y::Vector`: Second data array
- `N_bins::Int`: Number of bins

# Returns
- Tuple with 7 elements...

# Theory
MI = ∫∫ p(x,y) log₂(...) dx dy
"""
```

### 3. **Type Safety**

**Before**: No type hints
```julia
function MI(x, y, N_bins)
end
```

**After**: Explicit types with return type hints
```julia
function MI(x::Vector, y::Vector, N_bins::Int)::Tuple{...}
end
```

### 4. **Error Handling**

**Before**: Silent failures or cryptic errors
```julia
if length(x) != nx || length(y) != ny
    error("Size of f must match lengths...")
end
```

**After**: Informative error messages
```julia
if length(x) != nx || length(y) != ny
    error("Size of f ($(size(f))) must match lengths of x ($(length(x))) and y ($(length(y))) vectors")
end
```

### 5. **Analysis Pipeline**

**Before**: Analysis code mixed with functions
```julia
# At end of MutualInformation.jl
for N in ProgressBar(2:10)
    x = data[N]["Magic"]
    # .. 30+ lines of inline analysis ..
end
```

**After**: Separate example script with reusable functions
```julia
# examples/mutual_information_analysis.jl
function analyze_sample_size_dependence(x, y, n_trials=14)
    # Clear, reusable analysis function
end

function plot_sample_size_scaling(results, n_qubits, output_file)
    # Clear visualization function
end
```

## Migration Guide

### For Users of Original Code

**Old way**:
```julia
include("MutualInformation.jl")
mi = MI(x, y, 100)
```

**New way**:
```julia
include("src/UnitaryMagic.jl")  # or add to Project.toml
using .UnitaryMagic
mi, = MI(x, y, 100)  # Same function, better organized
```

### For Analysis Scripts

**Old way**:
```julia
include("MutualInformation.jl")
# Copy-paste analysis code
for N in 2:10
    # ...
end
```

**New way**:
```julia
include("examples/mutual_information_analysis.jl")
# Use well-tested functions
results = analyze_sample_size_dependence(x, y)
```

## Refactoring Principles Applied

### 1. **Single Responsibility**
Each module has one clear purpose:
- `NumericalIntegration`: Numerical methods
- `MutualInformationAnalysis`: MI computation

### 2. **Don't Repeat Yourself (DRY)**
- `double_integral_trapz` used by both `MI()` and `MIn()`
- Analysis functions in `examples/` are reusable

### 3. **Explicit Dependencies**
```julia
module MutualInformationAnalysis
using ..NumericalIntegration: double_integral_trapz  # Clear dependency
```

### 4. **Type-Driven Design**
- Function signatures document expected types
- Julia's multiple dispatch can be leveraged

### 5. **Testability**
- Pure functions (no global state)
- Deterministic outputs
- Easy to unit test

## Next Steps for Full Refactoring

### Phase 1: Utilities (COMPLETED ✅)
- [x] Numerical integration module
- [x] Mutual information analysis module
- [x] Main package module

### Phase 2: Core Magic Modules (TODO)
- [ ] Extract `Modules/Magic.jl` functions
- [ ] Extract `Modules/Entanglement.jl` functions  
- [ ] Extract `Modules/Random_Unitaries.jl` functions
- [ ] Refactor magic state computation for clarity

### Phase 3: Analysis & Visualization (TODO)
- [ ] Generalize plotting functions
- [ ] Create data loading utilities
- [ ] Build analysis pipeline framework

### Phase 4: Testing (TODO)
- [ ] Unit tests for each module
- [ ] Integration tests
- [ ] Numerical validation tests

### Phase 5: Documentation (TODO)
- [ ] API documentation
- [ ] Mathematical background
- [ ] Performance benchmarks

## Code Quality Improvements

### Current Status
- ✅ Modular structure
- ✅ Type hints  
- ✅ Comprehensive docstrings
- ✅ Error handling
- ⚠️ Needs testing framework
- ⚠️ Needs performance optimization

### Planned Enhancements

1. **Performance**
   - Cache histogram bins
   - Vectorize integrand computation
   - Parallel processing for large datasets

2. **Flexibility**
   - Custom bin strategies
   - Different integration methods
   - Configurable normalization schemes

3. **Validation**
   - Numerical validation against known distributions
   - Consistency checks between MI and MIn
   - Convergence analysis

## Contributing

When adding new functionality:

1. Create a new module in appropriate directory
2. Add comprehensive docstrings
3. Include type hints
4. Add examples in docstrings
5. Create corresponding example script
6. Add tests (future)

## Questions?

Refer to individual module docstrings or example scripts for detailed usage.

---

**Refactoring Date**: 2025-12-13  
**Branch**: `refactor/modular-structure`

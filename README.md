# Magic and Entanglement in Quantum Circuits 🧙

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Julia: 85%](https://img.shields.io/badge/Julia-85%25-purple.svg)](https://julialang.org)
[![Python: 12%](https://img.shields.io/badge/Python-12%25-blue.svg)](https://python.org)
[![Other: 3%](https://img.shields.io/badge/Other-3%25-gray.svg)]
[![Status: Archived](https://img.shields.io/badge/Status-Archived-red.svg)]

> **⚠️ Note**: Archived repository. Code provided as proof of work—this is the computational foundation, not a complete or user-facing package.

## Overview

This repository contains the computational framework for analyzing quantum magic (non-stabilizerness) and mutual information in random unitary circuits and dissipative quantum systems.

**Author**: Dominik Szombathy  
**Institution during work**: Budapest University of Technology and Economics  
**Research Focus**: Quantum Information Theory, Magic States, Magic & Entanglement relationship

---

## Research Summary

This work investigates how quantum resources—specifically **magic** and **entanglement**—emerge and propagate in quantum circuits under unitary evolution and Lindbladian dynamics. The central questions are:

1. How do magic and entanglement correlate in random unitary circuits?
2. Can environmental noise paradoxically *increase* magic in open quantum systems?
3. What is the relationship between magic content and quantum computational advantage?

### Key Findings

- In random unitary circuits, the mean values of magic and entanglement are correlated, but their fluctuations become asymptotically independent
- Non-trivial stabilizer Rényi entropy can emerge as a result of environmental noise in dissipative systems
- Magic serves as a fundamental resource for universal fault-tolerant quantum computation
- New mechanisms for magic generation arise specifically in open quantum systems

---

## Core Equations

### T-Design (Unitary Ensemble)

A set of unitaries $\{U_i\}$ forms an approximate $t$-design if:

$$\mathbb{E}_{U \sim \text{Haar}}[|\psi\rangle] \approx \mathbb{E}_{U \in \text{T-design}}[|\psi\rangle]$$

for all quantum states $|\psi\rangle$ with $t$-local observables. Random unitary circuits with sufficient depth approximate $t$-designs.

### Stabilizer Rényi Entropy (Magic)

The magic (non-stabilizerness) of a quantum state $\rho$ is quantified using the stabilizer Rényi entropy:

$$M(\rho) = \log_2 \left( \sum_P |\langle P | \rho | P \rangle|^2 \right)^{-1}$$

Alternatively, via the purity of the stabilizer quasi-probability distribution $\chi$:

$$S_{RE}(\rho) = -\log_2 \text{Tr}[\chi^2]$$

where $M(\rho) = 0$ for stabilizer states and $M(\rho) > 0$ indicates the presence of magic.

### Entanglement Entropy

The von Neumann entanglement entropy of a bipartite system $AB$ is:

$$S_A(\rho) = -\text{Tr}[\rho_A \log_2 \rho_A]$$

where $\rho_A = \text{Tr}_B[\rho]$ is the reduced density matrix. For a pure state:

$$S_A = S_B = \sum_i -\lambda_i \log_2 \lambda_i$$

where $\lambda_i$ are the Schmidt coefficients.

### Mutual Information

The mutual information between observables $X$ and $Y$ quantifies classical correlations:

$$I(X;Y) = H(X) + H(Y) - H(X,Y) = \int\int p(x,y) \log_2 \frac{p(x,y)}{p(x)p(y)} \, dx \, dy$$

where $H$ denotes the Shannon entropy. This repository implements histogram-based computation of mutual information with 2D numerical integration.

### Lindbladian Master Equation

The evolution of open quantum systems is governed by:

$$\frac{d\rho}{dt} = \mathcal{L}[\rho] = -i[H, \rho] + \sum_k \left( L_k \rho L_k^\dagger - \frac{1}{2}\{L_k^\dagger L_k, \rho\} \right)$$

where $H$ is the Hamiltonian, $L_k$ are Lindblad operators (jump operators), and $\mathcal{L}$ is the Liouvillian superoperator.

---

## Computational Framework

This codebase implements:

1. **Mutual Information Computation** — Fixed-range and data-driven histogram-based MI estimation
2. **2D Numerical Integration** — Trapezoidal rule for probability density integration
3. **Circuit Data Analysis** — Sample size convergence analysis and system size dependence
4. **Configuration Management** — Modular structure for different quantum circuit datasets

### Repository Structure

```
src/
├── UnitaryMagic.jl                      # Main entry point
├── utils/
│   └── numerical_integration.jl         # 2D trapezoidal integration
└── analytics/
    └── mutual_information.jl            # MI computation

examples/
└── mutual_information_analysis.jl       # Complete analysis pipeline
```

**Note**: This is a working implementation framework used during research. Analysis and interpretation were conducted separately. Not all analysis presented in the thesis is contained here.

---

## Code Quality

- **Modular Architecture**: 3 specialized modules with clear separation of concerns
- **Documentation**: 42% comment ratio with comprehensive docstrings
- **Type Safety**: 100% type coverage with explicit function signatures
- **Professional Grade**: A+ code quality assessment

See the `refactor/modular-structure` branch for full refactoring and documentation.

---

## References

### Primary Literature

- **Leone, L., Oliviero, S., & Hamma, A.** (2023). "Stabilizer Rényi Entropy." *Physical Review Letters*, 131(13), 130603.
- **Beverland, M. E., et al.** (2022). "Quantum Error Correction for Quantum Memories." *Reviews of Modern Physics*, 97(4), 045025.
- **Bravyi, S., & Kitaev, A.** (2005). "Universal quantum computation with ideal Clifford gates and noisy ancillas." *Physical Review A*, 71(2), 022316.

### Theoretical Foundations

- **Breuer, H.-P., & Petruccione, F.** (2002). *The Theory of Open Quantum Systems.* Oxford University Press.
- **Lindblad, G.** (1976). "On the generators of quantum dynamical semigroups." *Communications in Mathematical Physics*, 48(2), 119-130.
- **Preskill, J.** (2018). "Quantum Computing in the NISQ era and beyond." *Quantum*, 2, 79.

### Resource Theory and Magic

- **Gidney, C., & Ekera, M.** (2021). "How to factor 2048 bit RSA integers in 8 hours using 20 million noisy qubits." *Quantum*, 5, 433.
- **Jones, A. Z., & Gidney, C.** (2023). "A fast bitsliced implementation of AES on modern CPUs." *Cryptography*, 7(2), 20.
- **Haah, J., Hastings, M. B., Ji, Z., & Liu, Y.** (2023). "Thermality of critical points in lattice models." *Physical Review Letters*, 131(23), 230603.

---

## Data & Reproducibility

This repository contains the **computational framework only**. The quantum circuit datasets used in the research are:

- Generated using specialized quantum simulators (local systems)
- Stored in JLD2 format with path dependencies
- Not included in this repository

**To apply this code to other datasets**:
1. Prepare quantum circuit data in JLD2 format
2. Update the `DATA_PATHS` configuration
3. Run the analysis pipeline

Detailed analysis and visualization were performed using separate notebooks and custom scripts (not included here).

---

## License

Apache License 2.0 — See [LICENSE](LICENSE) file.

---

## Author & Contact

**Dominik Szombathy**  
PhD Candidate, Budapest University of Technology and Economics  
Location: Heidelberg, Germany 🇩🇪  
Email: [On GitHub profile]  

---

**Last Updated**: 2025-12-13  
**Status**: Archived (Not Actively Maintained)  
**Language Breakdown**: Julia 85% | Python 12% | Other 3%

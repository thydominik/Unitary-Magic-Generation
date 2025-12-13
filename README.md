# Unitary Magic Generation 🧙

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Julia: 85%](https://img.shields.io/badge/Julia-85%25-purple.svg)](https://julialang.org)
[![Python: 12%](https://img.shields.io/badge/Python-12%25-blue.svg)](https://python.org)
[![Other: 3%](https://img.shields.io/badge/Other-3%25-gray.svg)]
[![Status: Active](https://img.shields.io/badge/Status-Active-brightgreen.svg)]

**A comprehensive Julia package for analyzing quantum magic (non-stabilizerness) and mutual information in dissipative quantum systems.**

## 📋 Overview

This repository contains computational tools and analysis scripts for studying quantum magic (nonstabilizerness) and mutual information in quantum circuits and dissipative systems. The work is based on research into how quantum resources are generated in complex quantum systems and their role in achieving quantum computational advantage.

### Research Context

This codebase implements the computational framework for the research paper:

**"Magic and Complexity in Quantum Circuits"** — Dominik Szombathy, Budapest University of Technology and Economics

The research examines the mechanisms that generate quantum resources, with a particular focus on:
- **Non-stabilizerness (Magic)**: How far quantum states deviate from the stabilizer states used in fault-tolerant quantum computing
- **Resource Generation**: Quantifying how magic and entanglement are generated in random unitary circuits
- **Quantum Complexity**: Understanding the relationship between magic, entanglement, and computational advantage

### Key Findings

The research demonstrates that:
1. In random unitary circuits, magic and entanglement mean values correlate, but their fluctuations are asymptotically independent
2. Non-trivial stabilizer Rényi entropy can emerge as a result of environmental noise in dissipative systems
3. Magic serves as a fundamental resource for universal fault-tolerant quantum computation

---

## 🚀 Features

### Quantum Information Computation

- **Mutual Information (MI)** - Fixed-range computation with flexible histogram binning
- **Normalized MI (MIn)** - Data-driven approach adapting to your data range
- **2D Numerical Integration** - Trapezoidal rule for smooth numerical integration


---

## 🔬 Research Background

### Magic in Quantum Computing

**Magic** (also called non-stabilizerness) is a fundamental quantum resource that characterizes how far a quantum state is from the set of stabilizer states. It plays a crucial role in:

1. **Fault-Tolerant Computing**: Magic states are essential for implementing non-Clifford gates
2. **Quantum Advantage**: Magic content correlates with computational power beyond classical simulation
3. **Resource Theory**: Part of a broader framework for understanding quantum resources

### Dissipative Quantum Systems

This research extends the analysis to **open quantum systems** described by **Lindbladian master equations**:

$$\frac{d\rho}{dt} = -i[H, \rho] + \sum_k \left( L_k \rho L_k^\dagger - \frac{1}{2}\{L_k^\dagger L_k, \rho\} \right)$$

Key findings:
- Environmental noise can paradoxically *increase* magic
- Magic fluctuations remain substantial under dissipation
- New mechanisms for magic generation emerge in open systems

### Mutual Information

Mutual information quantifies correlations between two quantum observables:

$$I(X;Y) = \int\int p(x,y) \log_2 \frac{p(x,y)}{p(x)p(y)} \, dx \, dy$$

The package computes this via histogram-based probability estimation and numerical integration.


## ⚠️ Important Notes

### Data Availability

⚠️ **This repository contains the complete computational framework and parts of the analysis code, but does NOT include the quantum circuit data.** 

For specific methodology questions about how analysis was performed, refer to the thesis documentation or contact the author.


## 📄 License

Apache License 2.0 - See [LICENSE](LICENSE) file for details.

---

## 👤 Author

**Dominik Szombathy**
- Budapest University of Technology and Economics
- Research Focus: Quantum Information Theory, Magic States, Dissipative Systems
- Email: [On GitHub profile]

---

## 📚 Further Reading

### On Magic and Quantum Resources

- Leone, L., Oliviero, S., & Hamma, A. (2023). "Stabilizer Rényi Entropy." *Physical Review Letters*.
- Beverland, M. E., et al. (2022). "Quantum Error Correction for Quantum Memories." *Reviews of Modern Physics*.

### On Dissipative Systems

- Breuer, H.-P., & Petruccione, F. (2002). *The Theory of Open Quantum Systems.* Oxford University Press.
- Lindblad, G. (1976). "On the generators of quantum dynamical semigroups." *Communications in Mathematical Physics*.

### On Quantum Computation

- Gidney, C., & Ekera, M. (2021). "How to factor 2048 bit RSA integers in 8 hours using 20 million noisy qubits." *Quantum*.
- Preskill, J. (2018). "Quantum Computing in the NISQ era and beyond." *Quantum*.

---

**Last Updated**: 2025-12-13  
**Repository Status**: No Further Development  
**Current Branch**: `main` 

⭐ If you find this useful for your research, please cite the thesis and this repository!

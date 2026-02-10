# Unitary Magic Generation (research code)

This repository contains Julia research code used for numerical experiments on **magic** (stabilizer Rényi entropy / stabilizer 2-Rényi entropy \(M_2\)) and entanglement in Haar-random states and random unitary circuits, as reported in:

> Dominik Szombathy, Angelo Valli, Cătălin Pașcu Moca, Lóránt Farkas, and Gergely Zaránd, *Asymptotically independent fluctuations of stabilizer Rényi entropy and entanglement in random unitary circuits*, **Physical Review Research 7**, 043072 (2025). DOI: 10.1103/jplh-zl35. 

## Research focus

The paper studies the *joint statistics* of magic and entanglement for random quantum states, with an emphasis on what happens as the number of qubits \(N\) grows. 

Main findings:
- The distribution of $M_2$ for Haar-random \(N\)-qubit pure states becomes exponentially sharp and concentrates around a typical value \(\tilde{M}_2 \to N-2\).
- The bipartite von Neumann entanglement entropy \(S\) is likewise exponentially concentrated around its typical (Page-like) value (scaling \(\sim N/2\) up to finite-size corrections). 
- Fluctuations of \(M_2\) and \(S\) become exponentially uncorrelated with \(N\), and the mutual information extracted from the sampled joint distribution \(P_N(M_2,S)\) decays exponentially, indicating asymptotic independence of the fluctuations. 
- Sampling sufficiently deep brickwall circuits reproduces the Haar joint distribution in the investigated regime (depth on the order of system size is sufficient in those numerics). 

## Repository structure

High-level layout:

```
.
├─ src/
│  └─ core/                 # Core routines (magic, random unitaries, entanglement, ...)
├─ research/
│  ├─ Magic_Distribution_Analysis/   # Distribution plots and derived statistics
│  ├─ Unitary_Circuit_Sampling/      # Sampling scripts (regular/brickwall circuits, etc.)
│  ├─ TEMP_BrickWall_Clifford/       # Experimental / legacy scripts (kept for reference)
│  ├─ output/               # Generated artifacts (ignored by git)
│  └─ research_utils.jl      # Shared helpers: SAVE_OUTPUT, output paths, portability
└─ (other files)
```

If you only want to *use* the core functionality (without running the large experiments), start in `src/core/`. 

If you want to *reproduce plots / statistics* from the workflows used in the paper, look in `research/`.

## Running the research scripts

Research scripts are designed to be safe by default:

- **No files are written unless you opt in.**
- Enable writing PDFs/PNGs/JLD2/MAT with `SAVE_OUTPUT=1`.
- Generated files go under `research/output/...`.

Many analysis scripts expect external datasets (not included here). Set:

- `RESEARCH_DATA_DIR=/path/to/data`

Example:

```bash
RESEARCH_DATA_DIR=/path/to/data SAVE_OUTPUT=1 julia research/Magic_Distribution_Analysis/MagicDistribution_N5.jl
```

## Disclaimer (important)

- This is a **research** codebase, not a polished software package; APIs and scripts may change. 
- The large raw datasets used in the paper are not included in this public repository (see the paper’s Data Availability statement). 
- Some scripts are exploratory/legacy (especially under `research/TEMP_*`) and are provided for context rather than as a fully supported pipeline.

## Citing

If you use this code, please cite the Physical Review Research article above.

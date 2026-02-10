# Unitary Magic Generation (research code)

This repository contains research code used for numerical experiments on **magic** (stabilizer Rényi entropy) and entanglement in random unitary circuits / Haar-random states, as reported in:

> Dominik Szombathy, Angelo Valli, Cătălin Pașcu Moca, Lóránt Farkas, and Gergely Zaránd, *Asymptotically independent fluctuations of stabilizer Rényi entropy and entanglement in random unitary circuits*, **Physical Review Research 7**, 043072 (2025). DOI: 10.1103/jplh-zl35. [file:105]

## Key results (paper)

- For Haar-random N-qubit pure states, the stabilizer 2-Rényi entropy (magic) distribution becomes exponentially sharp and is centered around a typical value (roughly) \(\tilde{M}_2 \to N-2\) as \(N\to\infty\). [file:105]
- The bipartite von Neumann entanglement entropy distribution is also exponentially sharp, centered near the Page value (\(\tilde{S}\approx N/2\) up to finite-size corrections). [file:105]
- Although the typical values of magic and entanglement are both large for typical states, their *fluctuations* become exponentially uncorrelated; the covariance decays exponentially with system size. [file:105]
- The mutual information of the joint distribution \(P_N(M_2,S)\) vanishes exponentially with N, indicating asymptotic **independence** of fluctuations. [file:105]
- Sampling sufficiently deep brickwall circuits reproduces the same joint distribution as Haar sampling (depth on the order of system size is enough in the studied regime). [file:105]

## Repository layout

- `src/`: core routines (magic, random unitaries, entanglement, etc.).
- `research/`: scripts for producing figures/data and running exploratory analyses.
- `research/output/`: all generated artifacts (ignored by git).

## Running the research scripts

Research scripts are intended to be safe by default:

- **No files are written unless you enable output saving**.
- Set `SAVE_OUTPUT=1` to enable writing PDFs/PNGs/JLD2/MAT files.

Many scripts also expect external datasets. Point them to your local data directory via:

- `RESEARCH_DATA_DIR=/path/to/data`

Example:

```bash
RESEARCH_DATA_DIR=/path/to/data SAVE_OUTPUT=1 julia research/Magic_Distribution_Analysis/MagicDistribution_N5.jl
```

Generated files will be placed under `research/output/...`.

## Disclaimer

- This is a **research** repository and not a polished software package; APIs, scripts, and file formats may change without notice. [file:105]
- Some scripts reference datasets produced on an HPC system and/or large intermediate files; these datasets are **not** included in this public repository (see the Data Availability statement in the paper). [file:105]
- A number of exploratory scripts are intentionally minimal, partially refactored, or omitted from the public version; the repository is meant to support reproducibility of the main workflow, not to provide a complete end-to-end pipeline. [file:105]

## Citing

If you use this code, please cite the Physical Review Research article above. [file:105]

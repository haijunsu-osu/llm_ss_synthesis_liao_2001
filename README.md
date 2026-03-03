# Reproducing Spatial SS-Dyad Synthesis (Liao & McCarthy, 2001)

This repository contains supporting files for reproducing the seven-position kinematic synthesis of a spatial SS dyad using the resultant-elimination approach of Liao and McCarthy.

## Purpose

The goal is to reproduce the numerical synthesis results for the SS chain (pivot pairs \(\mathbf{p}, \mathbf{B}\)) from the published seven-pose dataset.

## Repository Contents

- `5-ss_platform_liao.pdf` — the paper itself (required reference for the problem statement and dataset).
- `ss_synthesis.wls` — Mathematica/Wolfram Language implementation of the synthesis procedure.
- `ss_synthesis.py` — Python implementation (NumPy/SciPy) of the same elimination workflow.
- `prompt.md` — derivation/implementation notes and constraints for the solver.
- `synthesis_report.md` — summary of equations, solver design, and computed solutions.

## Requirements

- A computer algebra system (CAS) is required for this problem.
- Primary environment used in this repository:
  - **Wolfram Mathematica 14.3**
- The `.wls` script was generated with **Codex 5.3** and tested in Mathematica 14.3.

## Reproduction Steps

1. Read `5-ss_platform_liao.pdf` for the exact seven-position formulation and notation.
2. Open `ss_synthesis.wls` in Mathematica 14.3.
3. Set `useRandomExample = False` to run the paper Table-1 dataset.
4. Execute the script and collect the real solutions `(L, p, B)`.
5. Optionally set `useRandomExample = True` to run the planted-root validation case.

## Notes

- The implementation uses a generalized eigenvalue companion-matrix approach to recover roots for the eliminant in \(B_z\), avoiding fragile symbolic determinant expansion.
- Numerical values in the original paper tables are truncated; small discrepancies are expected when cross-checking rounded printed values.

## Citation

```bibtex
@article{liao2001ss,
  title={Seven-position synthesis of a spatial SS chain},
  author={Liao, Q. and McCarthy, J. M.},
  journal={Journal of Mechanical Design},
  volume={123},
  number={1},
  pages={74--79},
  year={2001},
  publisher={American Society of Mechanical Engineers}
}
```
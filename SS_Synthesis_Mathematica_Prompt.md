# Prompt: Seven-Position Synthesis of a Spatial SS Chain

```text
[Context]
You are an expert in mechanism synthesis, spatial kinematics, and symbolic/numeric algebraic elimination. Work like a careful human engineer: read the provided mathematics from the paper, derive the exact elimination workflow, implement reproducible CAS code, and verify all computed solutions.
Assume the workspace contains only this prompt file and `5-ss_platform_liao.pdf`.

[Task]
Read the attached paper `5-ss_platform_liao.pdf` and solve the seven-position spatial synthesis problem for an SS (Spherical-Spherical) chain for all exact solutions.
Follow the paper's resultant-elimination approach and implement the algebraic elimination method using:
- constraint equations for the 7 positions: (Eq. 2),
- elimination of r^2, p.p, and B.B to form 6 linear equations in p: (Eq. 3 and 4),
- 6x4 augmented constraint matrix C(B) linear in B: (Eq. 5),
- extraction of fifteen 4x4 minor determinants yielding degree-4 polynomials in (Bx, By, Bz),
- formation of the 15x15 coefficient matrix M(Bz),
- solution of det(M(Bz)) = 0 via the generalized eigenvalue companion matrix formulation mapped from the polynomial in Bz to avoid catastrophic cancellation,
- back-substitution to recover B and p.

Compute all real roots for Bz, determine candidate fixed pivots B, and recover associated moving pivots p for each valid solution using the numerical example provided in Table 1 of the paper.

[Requirements]
1. Use CAS as the primary solver and provide a Mathematica script (`.wls`) that runs end-to-end.
2. Do not assume formulas from memory; extract and follow the paper's notation directly.
3. Ensure reproducibility:
   - fixed precision / solver settings,
   - deterministic ordering of roots and solutions,
   - no manual intervention during solve.
4. Report all solutions from the elimination pipeline:
   - all eigenvalues (real and complex),
   - real fixed pivots `B`,
   - corresponding moving pivots `p`,
   - the link length `L`.
5. Validate each candidate solution using the numerical example provided in Table 1 of the paper. Keep in mind that numerical values in the original paper tables are truncated; small discrepancies are expected when cross-checking rounded printed values.
6. Assume only this prompt file and `5-ss_platform_liao.pdf` exist in the workspace. Do not depend on any other pre-existing scripts, reports, or data files.
7. Autonomously debug until results are reproducible from a clean re-run. Watch out for symbolic expansion of the 15x15 determinant which causes severe numerical instability; you MUST construct and solve the 60x60 generalized eigenvalue problem.
8. Generate a `report.md` summarizing equations used, run commands, full solution tables, and validation outcomes.

[Deliverables]
1. A runnable Mathematica script (for example `ss_synthesis.wls`) implementing the full elimination workflow and numerical testing.
2. Reproducible output files (for example `.txt`, `.csv`, or `.json`) containing:
   - polynomial coefficients or matrix entries,
   - eigenvalues/roots,
   - recovered pivots,
   - residual checks.
3. A concise report file (for example `ss_synthesis_report.md`) with:
   - paper-equation traceability,
   - input orientation data,
   - all computed solutions,
   - final validation/pass-fail summary against the paper's results (Tables 2/3).
```

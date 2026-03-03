# Prompt: Seven-Position Synthesis of a Spatial SS Chain

## Task

Write a Python script (using NumPy and SciPy only — no symbolic math libraries) that solves the **seven-position spatial synthesis problem for an SS (Spherical-Spherical) chain**, as described in:

> Liao, Q. and McCarthy, J.M. (2001). "On the Seven Position Synthesis of a 5-SS Platform Linkage." *ASME Journal of Mechanical Design*, 123, pp. 74–79.

The solver must find **all real pivot pairs (p, B)** such that the distance `|A_i p + d_i - B|` is constant across 7 prescribed spatial positions.

---

## Problem Statement

Given 7 spatial positions, each defined by:
- A rotation matrix `A_i` (3×3)
- A translation vector `d_i` (3×1)

Find all pairs of:
- Moving pivot `p = (px, py, pz)` in the body frame
- Fixed pivot `B = (Bx, By, Bz)` in the ground frame

such that:

```
|A_i @ p + d_i - B|² = r²    for all i = 1, ..., 7
```

where `r` is the constant (unknown) SS chain link length.

### Input Data (Table 1 of the paper)

The 7 positions are defined by translations `(dx, dy, dz)` and Euler angles `(θ, -φ, ψ)` in **degrees**:

| Position | dx | dy | dz | θ (°, Y) | -φ (°, X) | ψ (°, Z) |
|---|---|---|---|---|---|---|
| 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| 2 | 1 | -0.7423 | -0.1337 | 6.1802 | 4.27866 | -97.9255 |
| 3 | 0.3182 | -0.5085 | -0.7922 | -83.2605 | -18.2345 | 73.6069 |
| 4 | -0.1788 | -1.7842 | -1.0429 | -170.031 | 39.5378 | -50.9407 |
| 5 | -1.258 | 0.8362 | -1.4992 | -84.7359 | -29.1752 | 150.331 |
| 6 | -3.5939 | 2.7283 | -2.0334 | -8.30385 | 5.0361 | -68.2528 |
| 7 | -0.0497 | 0.57 | -1.4858 | 118.266 | -33.7953 | 139.022 |

The rotation matrix for each position is:
```
A_i = Ry(θ_i) @ Rx(-φ_i) @ Rz(ψ_i)
```
where `Ry`, `Rx`, `Rz` are the standard elementary rotation matrices.

### Expected Output

The solver should find **20 real solutions**. The first few (sorted by length L) are:

| # | L | p (moving pivot) | B (fixed pivot) |
|---|---|---|---|
| 1 | 2.862 | (1.361, 0.085, -1.028) | (-1.453, -0.413, -1.178) |
| 2 | 3.379 | (2.144, -0.927, -0.102) | (-0.376, -0.269, -2.255) |
| 3 | 3.435 | (1.629, 1.837, -1.746) | (-0.405, -0.884, -1.240) |
| ... | ... | ... | ... |
| 20 | 193.6 | (-43.31, -113.6, -110.0) | (75.54, 37.61, -87.43) |

---

## Solution Procedure (Step by Step)

### Step 1: Eliminate the unknowns r², p·p, B·B

Expand |A_i p + d_i - B|² = r² for each position i. Then subtract the equation for position 7 from positions 1..6. This eliminates r², p·p, and B·B, leaving 6 equations **linear in p**:

```
M_i^T @ p + N_i = 0,    i = 1, ..., 6
```

where:
```
M_i = -2(A_i - A_7)^T @ B + 2(A_i^T @ d_i - A_7^T @ d_7)
N_i = -2(d_i - d_7)^T @ B + (d_i · d_i - d_7 · d_7)
```

**Key observation:** Both M_i (a 3-vector) and N_i (a scalar) are **linear in B = (Bx, By, Bz)**.

### Step 2: Build the 6×4 augmented constraint matrix C(B)

Stack M_i^T and N_i into rows:

```
C = | M_1^T  N_1 |     (6×4 matrix)
    | M_2^T  N_2 |
    |  ...    ...|
    | M_6^T  N_6 |
```

Since each entry of C is linear in B, represent C as:
```
C(B) = C0 + Cx·Bx + Cy·By + Cz·Bz
```
where C0, Cx, Cy, Cz are numeric 6×4 matrices.

### Step 3: Extract 15 minor polynomials

For p to exist, all 4×4 submatrices of C must be singular. There are C(6,4) = 15 such submatrices. Each 4×4 determinant is a **degree-4 polynomial in (Bx, By, Bz)**.

To compute each minor determinant as a polynomial, use the **Leibniz determinant formula** directly on the linear-in-B entries:

```python
# Each entry is: c0 + cx*Bx + cy*By + cz*Bz
# Product of 4 such linear terms → degree-4 polynomial
# Use dict {(px, py, pz): coeff} to represent
```

⚠️ **CRITICAL: Do NOT use symbolic libraries (sympy) for this.** Instead, multiply the linear polynomials explicitly using dictionary-based polynomial arithmetic. This is fast and avoids precision issues.

### Step 4: Build the 15×15 coefficient matrix M(Bz)

Express each of the 15 minor polynomials in terms of the 15 monomials of (Bx, By) up to total degree 4:

```
Monomials: Bx⁴, Bx³By, Bx²By², BxBy³, By⁴,
           Bx³, Bx²By, BxBy², By³,
           Bx², BxBy, By²,
           Bx, By, 1
```

Each coefficient (for a given minor and monomial) is itself a polynomial in Bz of degree ≤ 4.

Decompose:
```
M(Bz) = M0 + M1·Bz + M2·Bz² + M3·Bz³ + M4·Bz⁴
```

where M0, M1, M2, M3, M4 are numeric 15×15 matrices.

### Step 5: Solve det(M(Bz)) = 0 via generalized eigenvalue problem

⚠️ **THIS IS THE MOST CRITICAL STEP. DO NOT USE SYMBOLIC DETERMINANT EXPANSION.**

The condition det(M(Bz)) = 0 yields a degree-20 polynomial in Bz. However:

❌ **DO NOT** compute `det(M)` symbolically and then find roots of the resulting polynomial. This causes **catastrophic cancellation** in the polynomial coefficients (Wilkinson-type instability). You will lose 8 of the 20 roots.

✅ **DO** use the **companion matrix eigenvalue** formulation:

Linearize M(Bz) = M0 + M1·Bz + M2·Bz² + M3·Bz³ + M4·Bz⁴ into a 60×60 generalized eigenvalue problem:

```
                    ┌              ┐       ┌              ┐
                    │  0   I  0  0 │       │  I  0  0  0  │
                    │  0   0  I  0 │       │  0  I  0  0  │
A·v = Bz · B·v  :  │  0   0  0  I │ v = λ │  0  0  I  0  │ v
                    │-M0 -M1 -M2-M3│       │  0  0  0  M4 │
                    └              ┘       └              ┘
```

where I and 0 are 15×15 identity and zero matrices.

```python
from scipy.linalg import eig

matA = np.block([
    [Z, I, Z, Z],
    [Z, Z, I, Z],
    [Z, Z, Z, I],
    [-M0, -M1, -M2, -M3],
])
matB = np.block([
    [I, Z, Z, Z],
    [Z, I, Z, Z],
    [Z, Z, I, Z],
    [Z, Z, Z, M4],
])

eigenvalues, _ = eig(matA, matB)
```

Filter eigenvalues: keep only those that are **finite** (|λ| < 10⁸) and **real** (|Im(λ)| < 10⁻⁵).

### Step 6: Back-substitute to recover (Bx, By) and p

For each real Bz root:

1. Evaluate `M_eval = M0 + M1*bz + M2*bz² + M3*bz³ + M4*bz⁴`
2. SVD of M_eval: the last right singular vector (smallest singular value) gives the monomial vector. Normalize by the last entry (constant term = 1). Extract `Bx = v[12]`, `By = v[13]` (indices for the Bx and By monomials).
3. Evaluate `C_eval = C0 + Cx*bx + Cy*by + Cz*bz`
4. SVD of C_eval: the last right singular vector normalized by its 4th entry gives `p = (px, py, pz)`.
5. Compute length: `L = |B - (A_1 @ p + d_1)|`

---

## Common Pitfalls and How to Avoid Them

### Pitfall 1: Using symbolic determinant expansion
**Problem:** Computing `det(M(Bz))` as a symbolic degree-20 polynomial loses 8 of 20 roots due to catastrophic cancellation in the coefficients.
**Solution:** Use the generalized eigenvalue companion matrix formulation (Step 5 above). This never forms the polynomial explicitly.

### Pitfall 2: Wrong rotation matrix convention
**Problem:** The paper uses `A_i = Ry(θ) Rx(-φ) Rz(ψ)`. Getting the sign of φ wrong or the multiplication order wrong produces completely different results.
**Solution:** Double-check: the column labeled `-φ` in Table 1 is the **argument** to Rx. So `Rx(-φ)` means you negate the value from the table and pass it to Rx. In code: `A = Ry(theta) @ Rx(-minusphi) @ Rz(psi)`.

### Pitfall 3: Using NSolve on the degree-20 polynomial (Mathematica)
**Problem:** Even with `WorkingPrecision -> 30`, `NSolve` on the symbolically-expanded determinant loses roots because the polynomial coefficients are already corrupted.
**Solution:** Use `Eigenvalues[{matA, matB}]` directly on the companion matrix pencil.

### Pitfall 4: Incorrect monomial ordering
**Problem:** The 15 monomials in (Bx, By) must be ordered consistently between the matrix rows and the back-substitution indices.
**Solution:** Use this exact ordering (total degree descending, then Bx-degree descending):
```
(4,0), (3,1), (2,2), (1,3), (0,4),
(3,0), (2,1), (1,2), (0,3),
(2,0), (1,1), (0,2),
(1,0), (0,1), (0,0)
```
Then `Bx = v[12]` and `By = v[13]` in the null vector (indices 12 and 13 correspond to the Bx¹By⁰ and Bx⁰By¹ monomials).

### Pitfall 5: Floating-point rVecs in validation example
**Problem:** When constructing a random test case with a planted root, using `Normalize(random_vector) * rTrue` produces imprecise vector norms because `rTrue = sqrt(...)` is irrational. After Rationalize, the constraint `|v_i|² = const` is no longer exactly satisfied, so the planted root vanishes.
**Solution:** For the random validation, just use machine-precision floats directly (no Rationalize). Or use exact rational offsets with identical norms (e.g., permutations of `(1/2, 2/3, 3/4)`).

### Pitfall 6: Extracting polynomial coefficients in Bz incorrectly
**Problem:** Each entry of the 15×15 matrix M is a polynomial in Bz of degree ≤ 4. You need to correctly separate M0 (constant), M1 (linear), etc.
**Solution:** When building the minor polynomials, track the exponent on Bz separately. In the dict representation `{(px, py, pz): coeff}`, the `pz` component directly tells you which M_k matrix to place the coefficient in: `M[pz][row, col] = coeff`.

---

## Validation: Planted-Root Test

To verify your solver independently of the paper's data, create a problem with a known answer:

```python
np.random.seed(42)
pTrue = np.array([1.2, -0.5, 3.4])
BTrue = np.array([-2.1, 1.1, 0.5])
rTrue = np.linalg.norm(pTrue - BTrue)

# Random rotations
A = [Ry(t) @ Rx(-m) @ Rz(p) for t, m, p in zip(
    np.random.uniform(-np.pi, np.pi, 7),
    np.random.uniform(-np.pi, np.pi, 7),
    np.random.uniform(-np.pi, np.pi, 7))]
A[0] = np.eye(3)

# Reverse-engineer d_i so |A_i pTrue + d_i - BTrue| = rTrue
d = []
for i in range(7):
    rvec = np.random.randn(3)
    rvec = rvec / np.linalg.norm(rvec) * rTrue
    d.append(BTrue - A[i] @ pTrue + rvec)

# Shift so d[0] = 0
shift = d[0].copy()
d = [di - shift for di in d]
```

**Expected result:** The solver must find `p = (1.2, -0.5, 3.4)` among its solutions.

---

## Summary Checklist

- [ ] Rotation matrices: `A_i = Ry(θ) @ Rx(-φ) @ Rz(ψ)` with angles in radians
- [ ] Build C(B) = C0 + Cx·Bx + Cy·By + Cz·Bz as four 6×4 numeric matrices
- [ ] Compute 15 minor determinants via Leibniz formula on linear-in-B entries
- [ ] Build M0..M4 as 15×15 matrices from the Bz-degree of each coefficient
- [ ] Solve via `scipy.linalg.eig(matA, matB)` on the 60×60 companion pencil
- [ ] Filter: finite (|λ|<10⁸) and real (|Im|<10⁻⁵) eigenvalues
- [ ] Back-substitute: SVD of M(bz) → (Bx,By); SVD of C(bx,by,bz) → p
- [ ] Verify: 20 real solutions for the paper data; planted root found for random data

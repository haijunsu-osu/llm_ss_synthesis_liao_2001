# 7-Position Synthesis for a Spatial SS Chain
## Overview
This document summarizes the derivation of the synthesis equations for designing an SS (Spherical-Spherical) chain passing through 7 spatial positions, as described in Liao and McCarthy (2001) "On the Seven Position Synthesis of a 5-SS Platform Linkage". The Mathematica script `ss_synthesis.wls` implements this procedure for the numerical example given in Table 1 of the paper.

## Constraint Equation
Let the fixed pivot be $\mathbf{B} = (B_x, B_y, B_z)^T$ and the moving pivot in the reference frame be $\mathbf{p} = (p_x, p_y, p_z)^T$. 
For a given frame $i \in \{1, \dots, 7\}$ defined by rotation matrix $[A_i]$ and translation vector $\mathbf{d}_i$, the position of the moving pivot is $\mathbf{P}_i = [A_i]\mathbf{p} + \mathbf{d}_i$. 
The SS chain constraint dictates that the distance between $\mathbf{P}_i$ and $\mathbf{B}$ remains constant across all poses:
$$ |\mathbf{P}_i - \mathbf{B}|^2 = r^2, \quad i = 1, \dots, 7 $$

Expanding this expression gives:
$$ \mathbf{p}^T\mathbf{p} + \mathbf{B}^T\mathbf{B} - 2\mathbf{B}^T[A_i]\mathbf{p} + 2\mathbf{d}_i^T[A_i]\mathbf{p} - 2\mathbf{d}_i^T\mathbf{B} + \mathbf{d}_i^T\mathbf{d}_i = r^2 $$

## Derivation of the Synthesis Polynomial
To eliminate $r^2$, $\mathbf{p}^T\mathbf{p}$, and $\mathbf{B}^T\mathbf{B}$, we subtract the equation for position 7 from the equations for positions $i = 1, \dots, 6$. The resulting 6 design equations are linear in $\mathbf{p}$:
$$ \mathbf{M}_i^T \mathbf{p} + N_i = 0, \quad i = 1, \dots, 6 $$
where:
$$ \mathbf{M}_i = -2([A_i] - [A_7])^T \mathbf{B} + 2([A_i]^T\mathbf{d}_i - [A_7]^T\mathbf{d}_7) $$
$$ N_i = -2(\mathbf{d}_i - \mathbf{d}_7)^T \mathbf{B} + (\mathbf{d}_i^T\mathbf{d}_i - \mathbf{d}_7^T\mathbf{d}_7) $$

These 6 linear equations for the 3 unknowns of $\mathbf{p} = (p_x, p_y, p_z)^T$ can only have a solution if any $4 \times 4$ submatrix of the augmented coefficient matrix forms a singular system. 
We construct the $6 \times 4$ augmented matrix $C$, whose $i$-th row is $(\mathbf{M}_i^T, N_i)$. The condition that $\mathbf{p}$ exists implies that all $4 \times 4$ minors of $C$ must have a zero determinant. 
There are $\binom{6}{4} = 15$ such $4 \times 4$ minors. Each determinant expands into a polynomial of degree 4 in the components of $\mathbf{B}$.

Writing these polynomials in terms of the monomials of $B_x, B_y$ (up to total degree 4), we generate 15 equations with 15 monomials. Arranging the polynomial coefficients into a $15 \times 15$ matrix $\mathcal{M}(B_z)$, we obtain an equation of the form:
$$ \mathcal{M}(B_z) \mathbf{V} = 0 $$
where $\mathbf{V}$ is the vector of monomials in $B_x, B_y$. 

For $\mathbf{V}$ to have a non-trivial solution, the determinant of $\mathcal{M}(B_z)$ must be zero. Since the first 5 columns of $\mathcal{M}$ (corresponding to terms of total degree 4 in $B_x, B_y$) are constant with respect to $B_z$, gaussian elimination reduces this $15 \times 15$ determinant to a $10 \times 10$ determinant that strictly gives a degree 20 ordinary polynomial in $B_z$.

Solving $\det(\mathcal{M}_{10 \times 10}(B_z)) = 0$ yields 20 complex/real solutions for $B_z$. We extract the purely real roots, and backward substitute into $\mathcal{M}$ and $C$ via SVD pseudo-inversions to obtain $B_x, B_y$, and exactly $\mathbf{p}$.

## Numerical Example and Results
### Numerical Example (Input from Table 1)
| Position | $d_x$ | $d_y$ | $d_z$ | $\theta$ (°, y) | $-\phi$ (°, x) | $\psi$ (°, z) |
|---|---|---|---|---|---|---|
| 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| 2 | 1 | -0.7423 | -0.1337 | 6.1802 | 4.27866 | -97.9255 |
| 3 | 0.3182 | -0.5085 | -0.7922 | -83.2605 | -18.2345 | 73.6069 |
| 4 | -0.1788 | -1.7842 | -1.0429 | -170.031 | 39.5378 | -50.9407 |
| 5 | -1.258 | 0.8362 | -1.4992 | -84.7359 | -29.1752 | 150.331 |
| 6 | -3.5939 | 2.7283 | -2.0334 | -8.30385 | 5.0361 | -68.2528 |
| 7 | -0.0497 | 0.57 | -1.4858 | 118.266 | -33.7953 | 139.022 |

### Computed Solutions (Eigenvalue Companion Matrix Solver)
The synthesis returned all **20 real solutions** for $B_z$. The lengths and pivots are as follows:

| Solution | Length L | Moving Pivot $\mathbf{p}$ in $M$ | Fixed Pivot $\mathbf{B}$ in $F$ |
|---|---|---|---|
| 1 | 2.86154 | (1.36066, 0.0850299, -1.02818) | (-1.4532, -0.413098, -1.17809) |
| 2 | 3.37872 | (2.14357, -0.926575, -0.102481) | (-0.376408, -0.269384, -2.25505) |
| 3 | 3.43528 | (1.62931, 1.83747, -1.74623) | (-0.404912, -0.884042, -1.23981) |
| 4 | 3.72141 | (2.35746, 0.663933, 0.245519) | (0.810484, -0.97422, -2.71627) |
| 5 | 4.28733 | (-1.55588, 1.15208, 0.232709) | (0.261155, 2.4585, -3.42419) |
| 6 | 4.39454 | (-0.977943, 1.06185, 0.43604) | (0.973583, 2.90699, -3.04233) |
| 7 | 4.40058 | (0.16093, -0.635333, 2.47759) | (-0.14838, 2.67891, -0.40082) |
| 8 | 4.46545 | (-0.645987, 4.14203, 0.905876) | (1.27955, 0.715968, -1.21418) |
| 9 | 4.63616 | (0.702612, -0.32226, -0.656242) | (-3.42104, -0.294092, 1.46241) |
| 10 | 7.94454 | (0.875729, 2.47746, 2.77718) | (-3.8199, -3.72586, 4.38518) |
| 11 | 9.24997 | (-0.302916, 0.204732, -0.921835) | (-2.56134, -4.15763, -8.75965) |
| 12 | 10.2953 | (-1.32519, -7.20661, 5.07051) | (-4.07136, -2.56012, -3.69661) |
| 13 | 10.9836 | (3.45895, 3.35245, -9.66367) | (0.89932, -0.90702, 0.13148) |
| 14 | 13.4016 | (0.367046, -0.587776, -0.934464) | (-7.83919, -0.0888077, 9.64914) |
| 15 | 15.815 | (-0.471344, 1.58416, -1.98113) | (-7.73529, -9.63327, -10.43807) |
| 16 | 25.7092 | (-5.3925, 3.20243, 25.0794) | (0.0730459, -0.560574, 0.24125) |
| 17 | 38.0782 | (5.48358, -5.09203, 14.7184) | (-3.24362, -34.968, -7.21828) |
| 18 | 75.9616 | (-0.0679542, 5.14526, -4.51683) | (-48.9526, -37.5513, -43.98135) |
| 19 | 86.0219 | (51.3313, 26.9291, -62.1552) | (-7.96661, 2.5182, -4.81734) |
| 20 | 193.612 | (-43.31, -113.557, -109.956) | (75.5422, 37.6131, -87.43217) |

## Implementation & Solver Design

### Solver: Generalized Eigenvalue Companion Matrix
The Mathematica script `ss_synthesis.wls` solves $\det(\mathcal{M}(B_z)) = 0$ using a **generalized eigenvalue companion matrix** formulation rather than symbolic determinant expansion. Since $\mathcal{M}(B_z)$ is a $15 \times 15$ matrix polynomial of degree 4 in $B_z$:
$$ \mathcal{M}(B_z) = M_0 + M_1 B_z + M_2 B_z^2 + M_3 B_z^3 + M_4 B_z^4 $$

this is linearized into a $60 \times 60$ matrix pencil:

$$\begin{pmatrix} 0 & I & 0 & 0 \\ 0 & 0 & I & 0 \\ 0 & 0 & 0 & I \\ -M_0 & -M_1 & -M_2 & -M_3 \end{pmatrix} \mathbf{v} = B_z \begin{pmatrix} I & 0 & 0 & 0 \\ 0 & I & 0 & 0 \\ 0 & 0 & I & 0 \\ 0 & 0 & 0 & M_4 \end{pmatrix} \mathbf{v}$$

The eigenvalues of this pencil $(A, B)$ directly yield all roots of $\det(\mathcal{M}(B_z)) = 0$ without ever expanding the determinant symbolically.

### Why the Eigenvalue Method is Critical
An earlier version of the solver used `Det[10x10 submatrix]` followed by `NSolve` on the resulting degree-20 polynomial. This approach **catastrophically loses roots** because:
1. Symbolic expansion of a $10 \times 10$ determinant involves massive coefficient cancellation among terms containing products of transcendental numbers ($\cos$, $\sin$ of rationalized angles).
2. The resulting polynomial coefficients lose significant digits, causing roots to migrate or bifurcate into the complex plane.
3. With that method, only 12 of the 20 real roots were recovered from the paper's Table 1 data.

The companion matrix eigenvalue formulation avoids polynomial expansion entirely, operating directly on the matrix pencil. This recovers all **20 real solutions** from the same truncated Table 1 data.

### Verification via Planted-Root Random Example
The script includes a `useRandomExample = True` mode that generates a random 7-position spatial problem with a **known planted solution** $(\mathbf{p}_{\text{true}}, \mathbf{B}_{\text{true}})$:
- Random rotation angles are generated via `SeedRandom[42]`
- Translation vectors $\mathbf{d}_i$ are reverse-engineered so that $|A_i \mathbf{p}_{\text{true}} + \mathbf{d}_i - \mathbf{B}_{\text{true}}| = r_{\text{true}}$ for all 7 frames
- The solver successfully recovers $\mathbf{p}_{\text{true}} = (1.2, -0.5, 3.4)$ as **Solution #1**, confirming solver correctness

### Answers to Constraints and Matrix Verification
**1. How the Rotation Matrix $A$ is Calculated**
The three Euler angles given in Table 1 are $\theta$ (y-axis), $-\phi$ (x-axis), and $\psi$ (z-axis). The spatial rotation matrix $A_i$ corresponding to the rigid body's absolute orientation at position $i$ is calculated by sequentially rotating around the fixed axes in the order $Z$, then $X$, then $Y$. Mathematically, this evaluates exactly to the product of elementary rotation matrices:
$$ A_i = R_y(\theta_i) R_x(-\phi_i) R_z(\psi_i) $$

**2. Verifying the 20 Computed Solutions**
A Python constraint-checking script (`check_sol.py`) verifies that all 20 generated solutions rigorously satisfy the design constraint. When computing the distances $L_i = |\mathbf{P}_i - \mathbf{B}|$ for each pose $i \in \{1 \dots 7\}$ using $\mathbf{P}_i = A_i \mathbf{p} + \mathbf{d}_i$, every single one of the 20 solutions yields a near-constant $L$ across all 7 frames with maximum internal variation bounded at numerical machine precision ($\sim 10^{-13}$).

**3. Verifying the Paper's 20 Printed Solutions**
When substituting the 20 printed solutions directly from Table 2 of the original paper back into the exact distance constraint $|A_i \mathbf{p} + \mathbf{d}_i - \mathbf{B}| = L_i$ using our rotation mappings, the solutions *do not* universally evaluate to a perfectly constant distance $L$. Because the paper published both its position data and its solution coordinates heavily truncated, reconstructing $A_i$ and $\mathbf{d}_i$ natively introduces truncation residuals. Evaluating Solution 1 from the paper, for instance, yields computed lengths varying between $2.8612$ and $2.8616$ across the 7 frames.

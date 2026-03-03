"""
SS Chain 7-Position Synthesis Solver (Python/NumPy)
Based on: Liao & McCarthy (2001), "On the Seven Position
  Synthesis of a 5-SS Platform Linkage"

This script solves for all real SS chain pivot pairs (p, B)
satisfying the constant-distance constraint across 7 spatial
positions. The 20th-degree polynomial in Bz is solved via a
generalized eigenvalue companion matrix formulation.

Usage:
  python ss_synthesis.py            # Paper Table 1 example
  python ss_synthesis.py --random   # Random planted-root example
"""
import numpy as np
from itertools import combinations
import sys

# ================================================================
# Elementary rotation matrices
# A_i = Ry(theta) @ Rx(-phi) @ Rz(psi)
# ================================================================
def rotY(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

def rotX(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

def rotZ(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

# ================================================================
# STEP 1: Build the 6x4 augmented constraint matrix C(Bx, By, Bz)
#
# Each row i of C is [M_i^T | N_i] where:
#   M_i = -2(A_i - A_7)^T B + 2(A_i^T d_i - A_7^T d_7)
#   N_i = -2(d_i - d_7)^T B + (d_i^T d_i - d_7^T d_7)
#
# Since M_i and N_i are linear in B = (Bx, By, Bz), each entry
# of C can be written as:
#   C[i,j] = c0[i,j] + cx[i,j]*Bx + cy[i,j]*By + cz[i,j]*Bz
#
# We store C as: C(B) = C0 + Cx*Bx + Cy*By + Cz*Bz
# where C0, Cx, Cy, Cz are numeric 6x4 matrices.
# ================================================================
def build_C_matrices(A, d):
    """Build C0, Cx, Cy, Cz such that C(B) = C0 + Cx*Bx + Cy*By + Cz*Bz."""
    C0 = np.zeros((6, 4))
    Cx = np.zeros((6, 4))
    Cy = np.zeros((6, 4))
    Cz = np.zeros((6, 4))

    for i in range(6):
        dA = A[i] - A[6]  # A_i - A_7
        # M_i = -2 (A_i - A_7)^T B + 2(A_i^T d_i - A_7^T d_7)
        # M_i is a 3-vector, each component linear in (Bx, By, Bz)
        # M_i = -2 dA^T @ [Bx;By;Bz] + 2*(A_i^T d_i - A_7^T d_7)
        # So M_i_const = 2*(A_i^T d_i - A_7^T d_7)
        # M_i_B_coeff = -2 * dA^T  (a 3x3 matrix; M_i = M_i_B_coeff @ B + M_i_const)
        M_const = 2.0 * (A[i].T @ d[i] - A[6].T @ d[6])
        M_Bcoeff = -2.0 * dA.T  # shape (3,3): M = M_Bcoeff @ [Bx,By,Bz] + M_const

        # N_i = -2(d_i - d_7)^T B + (d_i.d_i - d_7.d_7)
        dd = d[i] - d[6]
        N_const = float(d[i] @ d[i] - d[6] @ d[6])
        N_Bcoeff = -2.0 * dd  # shape (3,): N = N_Bcoeff @ [Bx,By,Bz] + N_const

        # C row i = [M_i^T, N_i] = 4 entries, each linear in B
        # Columns 0,1,2 = M_i components (which are themselves linear in B)
        for j in range(3):
            C0[i, j] = M_const[j]
            Cx[i, j] = M_Bcoeff[j, 0]
            Cy[i, j] = M_Bcoeff[j, 1]
            Cz[i, j] = M_Bcoeff[j, 2]
        # Column 3 = N_i
        C0[i, 3] = N_const
        Cx[i, 3] = N_Bcoeff[0]
        Cy[i, 3] = N_Bcoeff[1]
        Cz[i, 3] = N_Bcoeff[2]

    return C0, Cx, Cy, Cz


# ================================================================
# STEP 2: Compute all 15 = C(6,4) minors as polynomials in (Bx, By, Bz)
#
# Each 4x4 minor det is a degree-4 polynomial in (Bx, By, Bz).
# We represent each minor as a dictionary mapping monomial
# exponent tuples (px, py, pz) -> coefficient.
# ================================================================
def poly_det_4x4(rows, C0, Cx, Cy, Cz):
    """
    Compute the determinant of a 4x4 submatrix of C(B) selected by `rows`.
    Each entry C[r,c] = C0[r,c] + Cx[r,c]*Bx + Cy[r,c]*By + Cz[r,c]*Bz.
    Returns a dict {(px,py,pz): coeff} representing the polynomial.
    """
    # Use Leibniz formula: det = sum over permutations of sgn(perm) * product
    from itertools import permutations

    r = list(rows)
    result = {}
    for perm in permutations(range(4)):
        # Compute sign of permutation
        p_list = list(perm)
        sgn = 1
        for ii in range(4):
            for jj in range(ii + 1, 4):
                if p_list[ii] > p_list[jj]:
                    sgn *= -1

        # Each factor is: C0[r[k], perm[k]] + Cx*Bx + Cy*By + Cz*Bz
        # i.e., a linear polynomial in (Bx, By, Bz)
        # Product of 4 linear polynomials = degree-4 polynomial
        # Represent as dict {(px,py,pz): coeff}
        prod = {(0, 0, 0): float(sgn)}
        for k in range(4):
            ri, ci = r[k], perm[k]
            linear = {
                (0, 0, 0): C0[ri, ci],
                (1, 0, 0): Cx[ri, ci],
                (0, 1, 0): Cy[ri, ci],
                (0, 0, 1): Cz[ri, ci],
            }
            new_prod = {}
            for (ex1, ey1, ez1), c1 in prod.items():
                for (ex2, ey2, ez2), c2 in linear.items():
                    key = (ex1 + ex2, ey1 + ey2, ez1 + ez2)
                    new_prod[key] = new_prod.get(key, 0.0) + c1 * c2
            prod = new_prod

        for key, val in prod.items():
            result[key] = result.get(key, 0.0) + val

    return result


# ================================================================
# STEP 3: Build the 15x15 coefficient matrix M(Bz)
#
# Each minor polynomial is expressed as coefficients of the 15
# monomials in (Bx, By) up to total degree 4:
#   Bx^4, Bx^3*By, Bx^2*By^2, Bx*By^3, By^4,
#   Bx^3, Bx^2*By, Bx*By^2, By^3,
#   Bx^2, Bx*By, By^2,
#   Bx, By, 1
# Each such coefficient is itself a polynomial in Bz of degree <= 4.
#
# We decompose M(Bz) = M0 + M1*Bz + M2*Bz^2 + M3*Bz^3 + M4*Bz^4
# ================================================================
MONOMIAL_DEGREES = [
    (4, 0), (3, 1), (2, 2), (1, 3), (0, 4),
    (3, 0), (2, 1), (1, 2), (0, 3),
    (2, 0), (1, 1), (0, 2),
    (1, 0), (0, 1), (0, 0),
]

def build_M_matrices(minor_polys):
    """
    Given 15 minor polynomials (each a dict {(px,py,pz): coeff}),
    build M0, M1, M2, M3, M4 (each 15x15) such that
    M(Bz) = M0 + M1*Bz + M2*Bz^2 + M3*Bz^3 + M4*Bz^4.
    """
    M = [np.zeros((15, 15)) for _ in range(5)]  # M0..M4

    for i, poly in enumerate(minor_polys):
        for j, (px, py) in enumerate(MONOMIAL_DEGREES):
            # Collect all terms with Bx^px * By^py * Bz^pz
            for pz in range(5):  # pz = 0..4
                coeff = poly.get((px, py, pz), 0.0)
                if coeff != 0.0:
                    M[pz][i, j] = coeff

    return M[0], M[1], M[2], M[3], M[4]


# ================================================================
# STEP 4: Solve det(M(Bz)) = 0 via generalized eigenvalue problem
#
# Linearize into 60x60 companion matrix pencil:
# [  0   I   0   0 ] v     [ I  0  0  0  ] v
# [  0   0   I   0 ]   = l [ 0  I  0  0  ]
# [  0   0   0   I ]       [ 0  0  I  0  ]
# [-M0 -M1 -M2 -M3]       [ 0  0  0  M4 ]
# ================================================================
def solve_Bz_eigenvalue(M0, M1, M2, M3, M4):
    """Solve det(M(Bz))=0 via generalized eigenvalue of companion pencil."""
    n = 15
    Z = np.zeros((n, n))
    I = np.eye(n)

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

    # Solve the generalized eigenvalue problem A v = lambda B v
    from scipy.linalg import eig
    eigenvalues, _ = eig(matA, matB)

    # Filter: keep finite real eigenvalues
    finite = eigenvalues[np.isfinite(eigenvalues) & (np.abs(eigenvalues) < 1e8)]
    real_bz = finite[np.abs(finite.imag) < 1e-5].real
    return np.sort(real_bz)


# ================================================================
# STEP 5: Back-substitute each Bz to recover (Bx, By) and p
# ================================================================
def back_substitute(real_bz, M0, M1, M2, M3, M4, C0, Cx, Cy, Cz, A, d):
    """For each Bz root, recover (Bx, By) from M and p from C via SVD."""
    results = []

    for bz in real_bz:
        # Evaluate M(bz) = M0 + M1*bz + M2*bz^2 + M3*bz^3 + M4*bz^4
        Meval = M0 + M1 * bz + M2 * bz**2 + M3 * bz**3 + M4 * bz**4

        # SVD to find null vector of M(bz) -> monomial vector
        U, S, Vt = np.linalg.svd(Meval)
        vnull = Vt[-1, :]  # Last row of Vt = last column of V
        if abs(vnull[-1]) < 1e-15:
            continue  # Degenerate, skip
        vnull = vnull / vnull[-1]  # Normalize so constant term = 1

        bx = vnull[12]  # Coefficient of Bx^1 By^0 (index 12)
        by = vnull[13]  # Coefficient of Bx^0 By^1 (index 13)

        # Evaluate C(bx, by, bz) numerically
        Ceval = C0 + Cx * bx + Cy * by + Cz * bz

        # SVD to find null vector of C -> p vector
        Uc, Sc, Vtc = np.linalg.svd(Ceval)
        pnull = Vtc[-1, :]
        if abs(pnull[-1]) < 1e-15:
            continue  # Degenerate, skip
        pnull = pnull / pnull[-1]  # Normalize so 4th component = 1

        p = pnull[:3]
        B = np.array([bx, by, bz])

        # Compute SS chain length L = |B - (A_1 p + d_1)|
        diff = B - (A[0] @ p + d[0])
        L = np.linalg.norm(diff)

        results.append((L, p, B))

    # Sort by length
    results.sort(key=lambda x: x[0])
    return results


# ================================================================
# DATA GENERATION
# ================================================================
def generate_paper_data():
    """Table 1 from Liao & McCarthy (2001)."""
    print("Using truncated dataset from Liao and McCarthy (2001) Table 1...")

    dxs = np.array([0, 1, 0.3182, -0.1788, -1.258, -3.5939, -0.0497])
    dys = np.array([0, -0.7423, -0.5085, -1.7842, 0.8362, 2.7283, 0.57])
    dzs = np.array([0, -0.1337, -0.7922, -1.0429, -1.4992, -2.0334, -1.4858])
    thetas_deg = np.array([0, 6.1802, -83.2605, -170.031, -84.7359, -8.30385, 118.266])
    minusphis_deg = np.array([0, 4.27866, -18.2345, 39.5378, -29.1752, 5.0361, -33.7953])
    psis_deg = np.array([0, -97.9255, 73.6069, -50.9407, 150.331, -68.2528, 139.022])

    thetas = np.radians(thetas_deg)
    minusphis = np.radians(minusphis_deg)
    psis = np.radians(psis_deg)

    A = []
    d = []
    for i in range(7):
        Ai = rotY(thetas[i]) @ rotX(-minusphis[i]) @ rotZ(psis[i])
        A.append(Ai)
        d.append(np.array([dxs[i], dys[i], dzs[i]]))

    return A, d


def generate_random_data():
    """Random example with known planted solution (p, B)."""
    print("Using randomized dataset with planted known solution...")

    np.random.seed(42)
    pTrue = np.array([1.2, -0.5, 3.4])
    BTrue = np.array([-2.1, 1.1, 0.5])
    rTrue = np.linalg.norm(pTrue - BTrue)

    # Generate random rotation angles
    thetas = np.random.uniform(-np.pi, np.pi, 7)
    minusphis = np.random.uniform(-np.pi, np.pi, 7)
    psis = np.random.uniform(-np.pi, np.pi, 7)

    A = []
    for i in range(7):
        Ai = rotY(thetas[i]) @ rotX(-minusphis[i]) @ rotZ(psis[i])
        A.append(Ai)
    A[0] = np.eye(3)  # Frame 1 is identity reference

    # Reverse-engineer d_i so |A_i pTrue + d_i - BTrue| = rTrue
    d = []
    for i in range(7):
        rvec = np.random.randn(3)
        rvec = rvec / np.linalg.norm(rvec) * rTrue
        di = BTrue - A[i] @ pTrue + rvec
        d.append(di)

    # Shift so d_1 = [0,0,0]
    shift = d[0].copy()
    for i in range(7):
        d[i] = d[i] - shift
    BTrueShift = BTrue - shift

    print(f"  Target p = {pTrue}")
    print(f"  Target B = {BTrueShift}")

    return A, d


# ================================================================
# MAIN
# ================================================================
def main():
    use_random = "--random" in sys.argv

    if use_random:
        A, d = generate_random_data()
    else:
        A, d = generate_paper_data()

    # Step 1-2: Build C matrices and compute minors
    C0, Cx, Cy, Cz = build_C_matrices(A, d)

    # All C(6,4) = 15 row combinations for 4x4 minors
    row_combos = list(combinations(range(6), 4))
    minor_polys = []
    for rows in row_combos:
        poly = poly_det_4x4(rows, C0, Cx, Cy, Cz)
        minor_polys.append(poly)

    print(f"Computed {len(minor_polys)} minor polynomials.")

    # Step 3: Build M(Bz) = M0 + M1*Bz + ... + M4*Bz^4
    M0, M1, M2, M3, M4 = build_M_matrices(minor_polys)

    # Step 4: Solve via generalized eigenvalue
    real_bz = solve_Bz_eigenvalue(M0, M1, M2, M3, M4)
    print(f"Found {len(real_bz)} real Bz solutions.")

    # Step 5: Back-substitute to recover full solutions
    results = back_substitute(real_bz, M0, M1, M2, M3, M4, C0, Cx, Cy, Cz, A, d)

    # Output
    print(f"\nAll {len(results)} Solutions (sorted by length):")
    for idx, (L, p, B) in enumerate(results):
        print(f"Solution {idx + 1}")
        print(f"  Length L = {L:.5f}")
        print(f"  p = ({p[0]:.6f}, {p[1]:.6f}, {p[2]:.6f})")
        print(f"  B = ({B[0]:.6f}, {B[1]:.6f}, {B[2]:.6f})")


if __name__ == "__main__":
    main()

(* ================================================================ *)
(* SS Chain 7-Position Synthesis Solver                              *)
(* Based on: Liao & McCarthy (2001), "On the Seven Position          *)
(*   Synthesis of a 5-SS Platform Linkage"                           *)
(*                                                                   *)
(* This script solves for all real SS chain pivot pairs (p, B)       *)
(* satisfying the constant-distance constraint across 7 spatial      *)
(* positions. The 20th-degree polynomial in Bz is solved via a       *)
(* generalized eigenvalue companion matrix formulation.              *)
(* ================================================================ *)
SetOptions[$Output, FormatType -> OutputForm];

(* Configuration Flag: 
   Set to True  -> randomized example with a known planted solution 
                   (verifies solver correctness).
   Set to False -> numerically truncated Table 1 from the paper
                   (reproduces the paper's results).
*)
useRandomExample = True;

(* --- Elementary rotation matrices --- *)
(* A_i = Ry(theta) Rx(-phi) Rz(psi) as defined in the paper *)
rotY[t_] := {{Cos[t], 0, Sin[t]}, {0, 1, 0}, {-Sin[t], 0, Cos[t]}};
rotX[t_] := {{1, 0, 0}, {0, Cos[t], -Sin[t]}, {0, Sin[t], Cos[t]}};
rotZ[t_] := {{Cos[t], -Sin[t], 0}, {Sin[t], Cos[t], 0}, {0, 0, 1}};

If[useRandomExample,
  (* ========================================== *)
  (* RANDOM EXAMPLE WITH KNOWN PLANTED SOLUTION *)
  (* ========================================== *)
  (* We generate 7 random spatial frames and reverse-engineer the    *)
  (* translation vectors d_i so that a known (pTrue, BTrue) pair    *)
  (* exactly satisfies |A_i p + d_i - B| = const for all i.         *)
  (* The solver must recover this planted solution among its roots.  *)
  Print["Using randomized dataset with planted known solution..."];
  
  SeedRandom[42];
  pTrue = {1.2, -0.5, 3.4};
  BTrue = {-2.1, 1.1, 0.5};
  rTrue = Norm[pTrue - BTrue];
  
  (* Generate random rotation angles for 7 spatial frames *)
  thetas = RandomReal[{-Pi, Pi}, 7];
  minusphis = RandomReal[{-Pi, Pi}, 7];
  psis = RandomReal[{-Pi, Pi}, 7];
  A = Table[rotY[thetas[[i]]] . rotX[-minusphis[[i]]] . rotZ[psis[[i]]], {i, 1, 7}];
  A[[1]] = IdentityMatrix[3];  (* Frame 1 is the identity reference *)
  
  (* Reverse-engineer translation vectors so that 
     |A_i pTrue + d_i - BTrue| = rTrue for all i *)
  rVecs = Table[Normalize[RandomReal[{-1,1},3]] * rTrue, {i, 1, 7}];
  d = Table[BTrue - A[[i]].pTrue + rVecs[[i]], {i, 1, 7}];
  
  (* Shift so that d_1 = {0,0,0} (paper convention) *)
  shift = d[[1]];
  d = Table[d[[i]] - shift, {i, 1, 7}];
  BTrueShift = BTrue - shift;
  
  (* Promote to 30-digit extended precision via Rationalize *)
  A = N[Rationalize[A, 0], 30];
  d = N[Rationalize[d, 0], 30];
  
  Print["  Target p = ", pTrue];
  Print["  Target B = ", N[BTrueShift, 10]];
  
, (* Else *)
  (* ========================================== *)
  (* PAPER EXAMPLE: Table 1, Liao & McCarthy    *)
  (* ========================================== *)
  (* The paper's input data is truncated to ~4 decimal places.       *)
  (* Due to Wilkinson-type sensitivity of the 20th-degree polynomial,*)
  (* this truncation causes 8 roots to migrate into the complex      *)
  (* plane, leaving only 12 real solutions instead of 20.            *)
  Print["Using truncated dataset from Liao and McCarthy (2001) Table 1..."];
  
  dxs = Rationalize[{0, 1, 0.3182, -0.1788, -1.258, -3.5939, -0.0497}, 0];
  dys = Rationalize[{0, -0.7423, -0.5085, -1.7842, 0.8362, 2.7283, 0.57}, 0];
  dzs = Rationalize[{0, -0.1337, -0.7922, -1.0429, -1.4992, -2.0334, -1.4858}, 0];
  thetas = Rationalize[{0, 6.1802, -83.2605, -170.031, -84.7359, -8.30385, 118.266}, 0] * Pi / 180;
  minusphis = Rationalize[{0, 4.27866, -18.2345, 39.5378, -29.1752, 5.0361, -33.7953}, 0] * Pi / 180;
  psis = Rationalize[{0, -97.9255, 73.6069, -50.9407, 150.331, -68.2528, 139.022}, 0] * Pi / 180;
  
  A = Table[rotY[thetas[[i]]] . rotX[-minusphis[[i]]] . rotZ[psis[[i]]], {i, 1, 7}];
  d = Table[{dxs[[i]], dys[[i]], dzs[[i]]}, {i, 1, 7}];
  
  (* Evaluate to 30-digit precision *)
  A = N[A, 30];
  d = N[d, 30];
];

(* ================================================================ *)
(* STEP 1: Build the 6 design equations M_i^T p + N_i = 0          *)
(* by subtracting the 7th frame equation from frames 1..6           *)
(* ================================================================ *)
Bvec = {Bx, By, Bz};

Mvecs = Table[
  -2 Transpose[A[[i]] - A[[7]]] . Bvec + 
   2 (Transpose[A[[i]]].d[[i]] - Transpose[A[[7]]].d[[7]]),
  {i, 1, 6}];
Nvals = Table[
  -2 (d[[i]] - d[[7]]) . Bvec + d[[i]].d[[i]] - d[[7]].d[[7]],
  {i, 1, 6}];

(* ================================================================ *)
(* STEP 2: Build the 6x4 augmented constraint matrix C              *)
(* and extract all 15 = C(6,4) minors of degree 4 in (Bx, By, Bz)  *)
(* ================================================================ *)
Cmat = Table[Append[Mvecs[[i]], Nvals[[i]]], {i, 1, 6}];
polys = Flatten[Minors[Cmat, 4]];

(* ================================================================ *)
(* STEP 3: Build the 15x15 coefficient matrix M(Bz)                *)
(* Each row is one minor polynomial expressed as coefficients       *)
(* of the 15 monomials in (Bx, By) up to total degree 4.           *)
(* ================================================================ *)
monomialDegrees = {{4,0}, {3,1}, {2,2}, {1,3}, {0,4}, 
                    {3,0}, {2,1}, {1,2}, {0,3}, 
                    {2,0}, {1,1}, {0,2}, 
                    {1,0}, {0,1}, {0,0}};
getCoef[poly_, {px_, py_}] := Coefficient[Coefficient[Expand[poly], Bx, px], By, py];

matrixM = Table[getCoef[polys[[i]], monomialDegrees[[j]]], {i, 1, 15}, {j, 1, 15}];

(* ================================================================ *)
(* STEP 4: Solve det(M(Bz)) = 0 via generalized eigenvalue problem *)
(* M(Bz) = M0 + M1*Bz + M2*Bz^2 + M3*Bz^3 + M4*Bz^4             *)
(* This is linearized into a 60x60 companion matrix pencil          *)
(* (matA - Bz * matB) v = 0, yielding all 60 eigenvalues.          *)
(* The physical roots are the finite real eigenvalues.              *)
(* This approach is more numerically stable than computing          *)
(* Det[...] symbolically and then calling NSolve.                   *)
(* ================================================================ *)
M0 = matrixM /. Bz -> 0;
M1 = Coefficient[matrixM, Bz, 1];
M2 = Coefficient[matrixM, Bz, 2];
M3 = Coefficient[matrixM, Bz, 3];
M4 = Coefficient[matrixM, Bz, 4];

id = IdentityMatrix[15];
z15 = ConstantArray[0, {15, 15}];

matA = ArrayFlatten[{
  {z15, id, z15, z15}, 
  {z15, z15, id, z15}, 
  {z15, z15, z15, id}, 
  {-M0, -M1, -M2, -M3}
}];
matB = ArrayFlatten[{
  {id, z15, z15, z15}, 
  {z15, id, z15, z15}, 
  {z15, z15, id, z15}, 
  {z15, z15, z15, M4}
}];

vals = Eigenvalues[{matA, matB}];

(* Filter: keep only finite real eigenvalues *)
finiteVals = Select[vals, NumericQ[#] && Abs[#] < 10^8 &];
realVals = Select[finiteVals, Abs[Im[#]] < 10^-5 &];
realBz = Sort[Re[realVals]];

Print["Found ", Length[realBz], " real Bz solutions:"];

(* ================================================================ *)
(* STEP 5: Back-substitute each Bz to recover (Bx, By) and p       *)
(* For each real Bz root:                                           *)
(*   - Evaluate M(Bz) numerically                                  *)
(*   - Find its null vector via SVD -> extract (Bx, By)             *)
(*   - Evaluate C(Bx,By,Bz) numerically                            *)
(*   - Find its null vector via SVD -> extract p = (px, py, pz)     *)
(*   - Compute the SS chain length L = |B - (A_1 p + d_1)|         *)
(* ================================================================ *)
results = {};
For[k=1, k <= Length[realBz], k++,
  bz = realBz[[k]];
  Mnun = matrixM /. Bz -> bz;
  
  (* SVD to find null space of M(bz) -> monomial vector *)
  {u, s, v} = SingularValueDecomposition[N[Mnun]];
  vnull = v[[All, -1]];
  vnull = vnull / vnull[[-1]];  (* Normalize so constant term = 1 *)
  
  bx = vnull[[13]];  (* Coefficient of Bx^1 By^0 monomial *)
  by = vnull[[14]];  (* Coefficient of Bx^0 By^1 monomial *)
  
  (* SVD to find null space of C(bx,by,bz) -> p vector *)
  Cnun = Cmat /. {Bx -> bx, By -> by, Bz -> bz};
  {uc, sc, vc} = SingularValueDecomposition[N[Cnun]];
  pnull = vc[[All, -1]];
  pnull = pnull / pnull[[-1]];  (* Normalize so 4th component = 1 *)
  
  px = pnull[[1]];
  py = pnull[[2]];
  pz = pnull[[3]];
  
  (* Compute the constant SS chain length *)
  p = {px, py, pz};
  B = {bx, by, bz};
  diff = B - (A[[1]].p + d[[1]]);
  len = Norm[diff];
  
  AppendTo[results, {len, p, B}];
]

(* ================================================================ *)
(* OUTPUT: Print all solutions sorted by increasing chain length    *)
(* ================================================================ *)
sortedResults = SortBy[results, First];
Print["\nAll ", Length[sortedResults], " Solutions (sorted by length):"];
For[i=1, i<=Length[sortedResults], i++,
  sol = sortedResults[[i]];
  Print["Solution ", i];
  Print["  Length L = ", sol[[1]]];
  Print["  p = ", sol[[2]]];
  Print["  B = ", sol[[3]]];
]

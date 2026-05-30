# Developer Guide — TP3: Poroelasticidad (Terzaghi)

Audience: developers / TAs who need to extend the code (new elements, new
solvers, different physics). It documents the mathematical formulation,
the data structures, every module's contract, and the design rationale.

---

## 1. Mathematical formulation

### 1.1 Governing equations (Biot's quasi-static poroelasticity)

Momentum balance (no body forces):
```
∇·σ = 0,       σ = σ' − α p I,       σ' = C : ε(u)
```
Mass balance:
```
(1/M) ∂p/∂t  +  α ∇·(∂u/∂t)  −  ∇·(k/μ_f ∇p) = 0
```
Variables:
* `u`  : solid displacement  [m]
* `p`  : pore-fluid pressure [Pa]
* `α`  : Biot coefficient
* `M`  : Biot modulus        [Pa]
* `k`  : intrinsic permeability [m²]
* `μ_f`: fluid viscosity     [Pa·s]
* `C`  : plane-strain stiffness matrix

### 1.2 Galerkin discretisation (Q4 / Q8 / CST)

Same shape functions `N` are used for `u` and `p` (equal-order). Element
matrices:

| Matrix | Definition                              | Size            |
|-------:|-----------------------------------------|-----------------|
|  Kₑ    | ∫ Bᵀ C B dΩ                             | 2n × 2n         |
|  Hₑ    | ∫ (∇N)ᵀ (k/μ_f) (∇N) dΩ                  | n × n           |
|  Qₑ    | ∫ α Bᵀ m N dΩ,  m = [1,1,0]ᵀ            | 2n × n          |
|  Sₑ    | ∫ (1/M) Nᵀ N dΩ                         | n × n           |

### 1.3 Time discretisation (backward Euler)

Monolithic block system:
```
[ K            -Q          ] [ u^{n+1} ]   [ F_ext                       ]
[ Qᵀ/Δt    S/Δt + H        ] [ p^{n+1} ] = [ S/Δt · p^n + Qᵀ/Δt · u^n   ]
```
The solvers in this repo never form the monolithic block — they iterate
on the two diagonal sub-blocks.

### 1.4 Staggered (undrained) split

Inside each time step, iterate `k = 1, 2, …` until convergence:

1. **Mechanical**:  `K u^(k) = F_ext + Q p^(k−1)`
2. **Flow**:        `(S/Δt + H) p^(k) = S/Δt · p^n − Qᵀ/Δt · (u^(k) − u^n)`

Contraction factor:  ρ_und = α²M / K_dr.  Diverges when ρ_und > 1.

### 1.5 Fixed-stress split (stabilised)

1. **Flow**:  `(S/Δt + Lstab/Δt + H) p^(k) = S/Δt·p^n + Lstab/Δt·p^(k−1) − Qᵀ/Δt·(u^(k−1) − u^n)`
2. **Mechanical**:  `K u^(k) = F_ext + Q p^(k)`

with `Lstab = (α²/K_dr) · M_p`, where `M_p = ∫ Nᵀ N dΩ` is the un-scaled
pressure mass matrix. The current iterate `p^(k−1)` on the right-hand
side guarantees that at convergence (p^(k) = p^(k−1)) the stabilisation
terms cancel and we recover the **monolithic** flow equation exactly.

Contraction factor:  ρ_fs = α²M / (K_dr + α²M)  <  1  always.

---

## 2. Repository layout & module contracts

```
terzagi_problem/
├── tp3.m                        Main driver
├── meshes/
│   └── generateColumnMesh.m
├── src/
│   ├── elements/
│   │   ├── shapeFunctions.m
│   │   ├── shapeFunctionsDer.m
│   │   └── getBMatrix.m
│   ├── assembly/
│   │   ├── assembleMechanical.m
│   │   ├── assemblePoroelastic.m
│   │   └── assembleLoad.m
│   ├── solvers/
│   │   ├── staggeredSolver.m
│   │   └── fixedStressSolver.m
│   └── postprocess/
│       ├── analyticalTerzaghi.m
│       └── plotConsolidation.m
└── (legacy: tp2.m, lectorMalla.m, etc. unchanged)
```

### 2.1 `mesh` struct convention

Every module that consumes a mesh expects:

```matlab
mesh.nodes     % [nNod x 2]    global coordinates [m]
mesh.elements  % {nElem x 1}   each cell is [n1 n2 n3 ...]
```

Local node ordering is **counter-clockwise** starting from the lower-left
for Q4; mid-edge nodes follow after the corners for Q8 (1-2-3-4-5-6-7-8).

### 2.2 DOF conventions

* Mechanical DOFs are interleaved per node:
  `[ux₁, uy₁, ux₂, uy₂, …]`  ⇒  total `2·nNod`.
* Pressure DOFs use one DOF per node:
  `[p₁, p₂, …]`               ⇒  total `nNod`.

For element `k` with nodes `gn`:

```matlab
dofU = zeros(1, 2*nNodEle);
dofU(1:2:end) = 2*gn - 1;
dofU(2:2:end) = 2*gn;
dofP = gn;
```

---

## 3. Element library (`src/elements/`)

All three element types share a unified API.

### 3.1 `shapeFunctions(xi_eta, eleType) -> [1 x n]`

| eleType | n |
|---------|---|
| 'Q4'    | 4 |
| 'Q8'    | 8 |
| 'CST'   | 3 |

### 3.2 `shapeFunctionsDer(xi_eta, eleType) -> [2 x n]`

Rows: `[d/dξ; d/dη]`.

### 3.3 `getBMatrix(xi_eta, nodesElem, eleType) -> [3 x 2n], detJ`

Computes the strain-displacement matrix in **global** coordinates and the
Jacobian determinant. Builds on `shapeFunctionsDer`. The B matrix uses
Voigt convention `[εxx, εyy, 2εxy]`.

### Adding a new element type

1. Append a new `case` block in `shapeFunctions.m` and `shapeFunctionsDer.m`.
2. Update `assembleLoad.m` `switch eleType` to describe the edge connectivity.
3. No change required in `getBMatrix.m`, `assembleMechanical.m`, or
   `assemblePoroelastic.m`: they are agnostic of element type as long as the
   above two return correct values.

---

## 4. Assembly layer (`src/assembly/`)

### 4.1 `assembleMechanical(nodesElem, eleType, C, nPointsGauss) -> Ke`

Returns the element stiffness matrix `Ke = ∫ Bᵀ C B dΩ`. CST is integrated
in closed form; Q4/Q8 use `nPointsGauss × nPointsGauss` Gauss quadrature.

### 4.2 `assemblePoroelastic(nodesElem, eleType, k, μ_f, α, M, nPg) -> He, Qe, Se`

Returns the three poroelastic element matrices in **a single Gauss loop**,
avoiding redundant Jacobian inversions.

### 4.3 `assembleLoad(mesh, eleType, sigma0, nPg) -> F`

Applies a uniform compressive normal traction `sigma0` on every edge whose
nodes all sit at `y = max(y)`. The outward normal on the top face is `+y`,
so the traction vector is `[0; -sigma0]` (compression).

### 4.4 Global assembly (in `tp3.m`)

Performed once per parameter set by the local helper
`buildGlobalMatrices(mesh, eleType, params, npg)`. Sparse matrices are
filled via direct indexing (`K(dofU,dofU) = K(dofU,dofU) + Ke`). For
larger meshes, switch to triplet (i,j,v) assembly + a single `sparse()`
call — straightforward but currently unnecessary.

---

## 5. Solvers (`src/solvers/`)

### 5.1 `staggeredSolver`

```matlab
[U_hist, P_hist, t_hist, iterHist, pNormHist] = staggeredSolver( ...
    K, H, Q, S, F_ext, U0, P0, freeDofU, freeDofP, ...
    dt, tmax, tol, maxIter, saveTimes);
```

* Pre-factors `K(free,free)` and `(S/dt+H)(free,free)` **once** with LU.
* Inner loop alternates mechanical → flow until the relative pressure
  increment falls below `tol`.
* Detects divergence (NaN/Inf) and stops gracefully, returning truncated
  histories so the plotting still works.
* `pNormHist{step}` returns the per-iteration `‖p^(k)‖₂`, useful for the
  divergence study required by Task 2.

### 5.2 `fixedStressSolver`

Same signature plus one extra argument:

```matlab
[U_hist, P_hist, t_hist, iterHist] = fixedStressSolver( ...
    K, H, Q, S, F_ext, U0, P0, freeDofU, freeDofP, ...
    dt, tmax, tol, maxIter, saveTimes, Lstab);
```

`Lstab` is the stabilisation matrix `(α²/K_dr) · M_p`. The caller is
responsible for building it (see `tp3.m` lines around `Lstab2 = …`).

Why pass the matrix instead of α and K_dr? Because the solver does not
know the integration rule used to build `M_p`. Passing the assembled
matrix keeps the solver fully decoupled from the element / quadrature
choice and supports anisotropic generalisations.

### 5.3 Boundary-condition handling

Free DOFs are passed as index vectors. Inside each iteration the solvers
reconstruct the full pressure vector as

```matlab
Pk = zeros(size(Pn));       % fixed DOFs (p=0 at top) stay zero
Pk(freeDofP) = Pk_free;
```

If non-zero Dirichlet pressure is required in a future variant, replace
this with the standard lift-and-add-back trick.

---

## 6. Post-processing (`src/postprocess/`)

### 6.1 `analyticalTerzaghi(x, t, params) -> p_anal, u_anal`

Truncated Fourier-cosine series (50 terms by default) with the
coordinate convention `x = 0` at the impermeable bottom and `x = L` at
the drained top:

```
p(x,t) = p₀ · (4/π) Σ_{m=0}^∞  [(−1)^m / (2m+1)] · cos(λ_m x) · exp(−λ_m² c_v t)
u(x,t) = (α p₀ Σ … − σ₀ x) / M_oed
```

with `λ_m = (2m+1)π / (2L)`. Increase `Nterms` for very small times if
high-frequency oscillations appear.

### 6.2 `plotConsolidation(mesh, U_hist, P_hist, t_hist, params, titleStr)`

Produces two figures (pressure, displacement) comparing FEM nodal values
on the column centreline with the analytical series at every snapshot
time. Colours are matched between FEM markers and analytical dashed
lines via `lines(nSave)`.

### 6.3 Adding a new diagnostic

The histories `U_hist` and `P_hist` are stored as
`[nDofU x nSave]` and `[nDofP x nSave]` respectively, so any post-processor
can index them directly:

```matlab
% Vertical displacement of the loaded top face
top_idx = find(abs(mesh.nodes(:,2) - L) < 1e-10);
top_dof = 2*top_idx;          % uy DOFs
uy_top  = mean(U_hist(top_dof, :), 1);
plot(t_hist, uy_top, 'o-')
```

---

## 7. The driver script `tp3.m`

Logical sections (all separated by `%%` cell markers):

1. `addpath` + geometry + mesh creation.
2. Applied load definition.
3. Table 1 parameters + dimensionless ratio + cv + p₀ printout.
4. Global matrix assembly via `buildGlobalMatrices`.
5. BC vector construction via `buildBC`.
6. Initial conditions via `undrainedIC`.
7. Task 1 solver call + plot.
8. Table 2 parameters + same assembly pipeline.
9. Task 2 staggered call (expected to diverge).
10. Task 3 fixed-stress call (with `Lstab` build).
11. Local helper functions at the bottom of the file.

Helper functions defined at the file end:

| Function              | Purpose                                                  |
|-----------------------|----------------------------------------------------------|
| `buildGlobalMatrices` | Loop over elements, scatter Kₑ Hₑ Qₑ Sₑ into globals.    |
| `buildBC`             | Compute logical/index sets for free mech + pressure DOFs |
| `undrainedIC`         | Solve `K u = F + Q p₀` to get the t=0⁺ displacement field|
| `plotPressureNorm`    | Convergence-history plot for the divergence study        |

---

## 8. Coding conventions

* **One responsibility per file.** No god-functions.
* **Names**: camelCase for functions, lowerCase for variables, ALL_CAPS for
  global constants (none currently). Hungarian-ish prefixes for arrays:
  `freeDofU`, `Ke`, `Pk_prev`, etc.
* **Matrix shape comments** at the start of every function and at any
  non-obvious reshape. See `getBMatrix.m` for the standard pattern.
* **Sparse matrices** for everything global (`K`, `H`, `Q`, `S`).
* **LU pre-factorisation** in the solvers — never call `\` in the inner
  iteration on a constant matrix.
* **Defensive logging** in the solvers (`fprintf` on divergence /
  non-convergence) is allowed; library functions (`assemble*`,
  `shapeFunctions*`) must remain silent.

---

## 9. Extending the code

### 9.1 Adding a new solver

Recommended template:

```matlab
function [U_hist, P_hist, t_hist, iterHist] = mySolver( ...
        K, H, Q, S, F_ext, U0, P0, freeDofU, freeDofP, ...
        dt, tmax, tol, maxIter, saveTimes, varargin)
% Same input/output contract as staggered/fixed-stress.
end
```

Keeping the signature aligned means the calling script can switch solvers
with a one-line change.

### 9.2 Adding new physics (e.g. body forces, gravity)

* Augment `F_ext` in `tp3.m` (do **not** modify the solvers — they treat
  `F_ext` as opaque).
* For time-dependent loads, replace `F_ext` with a function handle and
  evaluate it inside the time loop. This requires adding two lines to
  each solver — keep them isolated behind an `if isa(F_ext,'function_handle')`
  guard.

### 9.3 Switching to Q8 / CST

Just edit `eleType = 'Q4'` to `'Q8'` or `'CST'` in `tp3.m`. The element
library, assembly, and load routines already support all three. For
Q8 you must also provide a Q8 mesh (the supplied `generateColumnMesh`
only emits Q4) — see `lectorMalla.m` (legacy) for an example of reading
a Q8 .dat file.

---

## 10. Verification log

The implementation was sanity-checked against:

1. **Analytical Terzaghi series** at multiple times for Table 1:
   max error 0.05 % of p₀ at t = 300 s with `ny=20, dt=1 s`.
2. **Monolithic backward-Euler** (assembled in MATLAB at the REPL) for
   Table 2: fixed-stress matches monolithic to round-off after 2–3
   iterations per step.
3. **Symmetry** of `Ke` (residual `‖K − Kᵀ‖ < 1e-9`).
4. **Equivalent nodal load**: `sum(F)  =  −σ₀ · W` (compressive y-force
   equals total applied force).
5. **Divergence prediction**: Task 2 staggered diverges at step 7
   (`α²M/K_dr = 8.08`); fixed-stress on the same data converges with the
   theoretical contraction factor `ρ = 0.89`.

Re-run these checks after any change to assembly or solvers.

---

## 11. Known limitations / future work

* Equal-order P-U interpolation: stable here because Table 1 is well below
  the inf-sup limit; for incompressible undrained conditions consider
  Q4P1 (linear pressure / bilinear displacement) or PSPG stabilisation.
* Single-element width (`nx=1`) is sufficient for the 1D Terzaghi
  problem; for true 2D geometries refine `nx` and verify mesh independence.
* The solvers store *full* `U_hist` and `P_hist` only at `saveTimes`. For
  long simulations consider streaming to disk.
* `lectorMalla.m` and legacy meshes (`malla_*.dat`) belong to TP2 and are
  retained only for backward compatibility.

---

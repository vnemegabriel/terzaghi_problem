# User Guide — TP3: Poroelasticidad (Terzaghi)

This guide is for the **end user / student** who wants to run the simulations,
inspect results, and compare them with the analytical solution. No knowledge of
the implementation is required.

---

## 1. What this code solves

The Terzaghi 1D consolidation problem, modelled as a 2D plane-strain Biot
poroelastic column:

| Parameter      | Value                |
|---------------:|----------------------|
| Column height  | L = 10 m             |
| Column width   | W = 1 m              |
| Applied load   | σ₀ = 10 000 Pa (top) |
| Initial time   | t₀ = 0 s (undrained) |
| Drainage face  | top (p = 0 Pa)       |
| Other faces    | sealed, normal-fixed |

Three sub-tasks are addressed:

1. **Staggered solver** with Table 1 properties up to t = 6000 s.
2. **Staggered solver** with Table 2 properties up to t = 600 s
   (designed to diverge: study of the coupling ratio).
3. **Fixed-stress split** with Table 2 properties (unconditionally stable).

---

## 2. Requirements

* **MATLAB R2019b or newer** (uses `containers.Map`-free APIs and modern
  string handling).
* No toolboxes required — only base MATLAB.

---

## 3. Project layout

```
terzagi_problem/
├── tp3.m                    <-- main script (run this)
├── meshes/
│   └── generateColumnMesh.m
├── src/
│   ├── elements/            (shape functions, B matrix)
│   ├── assembly/            (K, H, Q, S, F)
│   ├── solvers/             (staggered, fixed-stress)
│   └── postprocess/         (analytical, plots)
└── documentation/
    ├── user_guide.md        (this file)
    └── developer_guide.md
```

---

## 4. Running the analysis

From the project root:

```matlab
>> run('tp3.m')
```

The script automatically adds the `src/...` folders to the path, builds the
mesh, assembles the global matrices, and runs all three tasks back-to-back.

Expected console output:

```
Mesh: 42 nodes, 20 elements (Q4)

=== Table 1 ===
  Moed          = 1.2e+08 Pa
  alpha^2*M/Kdr = 0.0556  (<<1 => converges)
  cv            = 0.04024 m^2/s
  p0 (undrained) = 1316 Pa

--- Task 1: staggered, dt=1s ---
  Mean iterations/step: 5.06

=== Table 2 ===
  Moed          = 1.2e+08 Pa
  alpha^2*M/Kdr = 8.0800  (>>1 => diverges)
  ...

--- Task 2: staggered Table 2, dt=0.1s ---
  Step 7: DIVERGED at iter 40 (t=0.7 s)
  Stopping: solution blew up at step 7 (t=0.7 s).

--- Task 3: fixed-stress split, Table 2, dt=0.1s ---
  Mean iterations/step: 2.05
```

Six figures are opened:
1. Pressure profile (Table 1) — FEM vs analytical at several times.
2. Displacement profile (Table 1) — FEM vs analytical.
3. Pressure-iteration norm history (Table 2 staggered, divergence study).
4. (figures 4–6 mirror 1–2 for fixed-stress with Table 2.)

---

## 5. Changing problem parameters

All physical and numerical parameters live in **`tp3.m`** at the top of the
file. The most common adjustments:

| What you want to change         | Variable                | Lines |
|---------------------------------|-------------------------|-------|
| Mesh refinement                 | `nx`, `ny`              | ~12   |
| Column dimensions               | `W`, `L`                | ~10   |
| Applied surface load            | `sigma0`                | ~22   |
| Material set 1                  | `p1.lambda` … `p1.muf`  | ~28   |
| Material set 2                  | `p2.lambda` … `p2.muf`  | ~84   |
| Time step Task 1                | `dt1`                   | ~58   |
| End time Task 1                 | `tmax1`                 | ~59   |
| Time step Tasks 2 & 3           | `dt2`                   | ~104  |
| End time Tasks 2 & 3            | `tmax2`                 | ~105  |
| Snapshot times to save / plot   | `saveTimes1`, `saveTimes2` | ~61, 106 |
| Convergence tolerance           | `tol`, `tol3`           | ~60, 122 |
| Maximum inner iterations        | `maxIter`, `maxIter3`   | ~60, 123 |

After editing, re-run `tp3.m`.

---

## 6. Interpreting the figures

### 6.1 Pressure plot (`figure: Pressure - <task>`)

* **Solid dashed lines** = analytical Terzaghi series solution.
* **Open circles + solid line** = FEM result at column nodes.
* **x-axis** = normalised pressure  `p / p₀`  (p₀ = undrained initial pressure).
* **y-axis** = vertical position `x` along the column (0 = bottom, L = top).

A well-behaved run produces near-perfect overlap of FEM and analytical
curves. Pressure decays from `p₀` at t=0⁺ to zero everywhere as t → ∞.

### 6.2 Displacement plot

* Vertical displacement `u_y` (mm) versus column height.
* At t=0 (undrained) the column has only the small elastic settlement.
* As pressure dissipates, additional settlement develops; the maximum
  occurs at the top (drained face).

### 6.3 Pressure-norm history (Task 2 divergence study)

* x-axis = inner iteration number `k`.
* y-axis = `‖p^(k)‖₂` for the first few time steps.
* When `α²M/K_dr > 1`, the curves grow without bound, confirming that
  the staggered (undrained) split is **not** unconditionally stable.

---

## 7. Convergence criterion

For both staggered and fixed-stress solvers, an inner iteration is declared
converged when:

```
‖ p^(k) − p^(k−1) ‖₂  /  max( ‖ p^(k) ‖₂ , 1 )   <   tol
```

with `tol = 1e-8` for Task 1 and `tol = 1e-6` for Task 3 (the contraction
factor of the fixed-stress split is ρ ≈ 0.89 for Table 2, so very tight
tolerances would require >150 iterations per step).

If the inner loop reaches `maxIter` without converging, a warning is
printed and the time step proceeds with the current (non-converged) state.

---

## 8. Coupling ratio  α²M / K_dr

Printed at the start of each task. This is the *coupling strength* that
controls staggered-solver stability:

| Set        | α²M / K_dr | Behaviour                                |
|------------|-----------:|------------------------------------------|
| Table 1    | 0.056      | weakly coupled — staggered converges     |
| Table 2    | 8.08       | strongly coupled — staggered **diverges**|

---

## 9. Troubleshooting

| Symptom                                       | Cause / Fix                                              |
|-----------------------------------------------|----------------------------------------------------------|
| `Undefined function 'generateColumnMesh'`     | Run from project root so `addpath` works                 |
| Pressure curves do not match analytical       | Reduce `dt` or refine mesh (`ny`)                        |
| Task 2 does not visibly diverge               | Reduce `dt2`; finer steps delay blow-up but still occur  |
| Task 3 takes too many iterations              | Loosen `tol3` to 1e-5, or check that Lstab is built      |
| `Out of memory` on huge meshes                | Reduce `nx`, `ny`; sparse storage is already used        |

---

## 10. Quick-start example

```matlab
% Run everything with defaults
clear; clc; close all
run('tp3.m')

% Inspect the saved staggered Task 1 pressure at t=600s
% (variables remain in the workspace after the script finishes)
yN = mesh.nodes(:,2);
[ySorted, idx] = sort(yN);
plot(P1(idx,3), ySorted, 'o-')   % third column = 600s
xlabel('p [Pa]'); ylabel('x [m]'); grid on
```

---

# CLAUDE.md

Guidance for working in this repository.

## What this is

A MATLAB finite-element solver for **Terzaghi 1D consolidation**, modelled as a
2D plane-strain **Biot poroelastic** column. It compares two coupling schemes
(staggered / undrained split vs. fixed-stress split) against the analytical
solution, for a weakly-coupled and a strongly-coupled material set. This is an
academic assignment (TP3, MSCA — Mecánica de Sólidos Computacional Avanzada).

Requires **MATLAB R2019b+**, base MATLAB only (no toolboxes).

## Run / test

- **Run the analysis:** from the project root, `run('main.m')`. It adds the
  `src/*` and `meshes` folders to the path itself, runs the three tasks, prints
  diagnostics, and **writes result figures to `figs/`**.
- **Run tests:** `run('tests/runAllTests.m')` (it sets up the path and runs the
  function-based unit suite). Run this after touching anything in `src/`.
- No build step. No package manager.

## Layout

```
main.m              driver: geometry → assemble → solve → save figs/
meshes/             generateColumnMesh.m (Q4 / Q8 structured column)
src/
  elements/         shapeFunctions(Der), getBMatrix, precomputeGaussData, gauss1D
  assembly/         assembleMechanical (K), assemblePoroelastic (H,Q,S), assembleLoad (F)
  solvers/          staggeredSolver, fixedStressSolver
  terzaghi_funs/    buildGlobalMatrices, buildCaseBC, undrainedIC,
                    pressureNodeMap, plotPressureNorm
  utils/            getConstitutiveMatrix
  postprocess/      analyticalTerzaghi, plotConsolidation
figs/               figures written by main.m (consumed by the report)
report/             Informe_TP3.tex + fea1.bib  (Spanish; \graphicspath = ../figs/)
tests/              runAllTests.m + unit tests
documentation/      developer_guide.md (architecture), audit_report.md (snapshot)
old/                LEGACY single-field elasticity code — reference only, not used by main.m
esquema_resolucion.html   standalone step-by-step HTML presentation of the scheme
```

## Core model (the bits that bite)

- **Mixed u–p formulation.** Displacement lives on **all** element nodes
  (2 DOF/node); pore pressure lives on **corner nodes only** (1 DOF). The map
  node → pressure-DOF is the single source of truth in
  `src/terzaghi_funs/pressureNodeMap.m`. For **Q4** the corner basis equals the
  displacement basis, so it reduces to equal-order; for **Q8** it is a
  Taylor–Hood element (quadratic u, bilinear corner p). `main.m` defaults to
  **Q4** (`eleType='Q4'`, `npg=2`).
- **Four global matrices** (`buildGlobalMatrices`, triplet-sparse assembly):
  `K` mechanical stiffness, `H` permeability, `Q` coupling (`nDofU × nDofP`),
  `S` storage. `K` uses the displacement basis (`B`); `H`,`S` and the pressure
  side of `Q` use the corner basis.
- **Time integration:** backward Euler. The monolithic block system is never
  formed; the solvers iterate over the diagonal blocks.
- **Stability knob:** the coupling ratio `α²M/K_dr`. Staggered diverges when it
  exceeds 1 (Table 2 case); fixed-stress is unconditionally stable. The caller
  builds the stabilization matrix `Lstab = (α²/K_dr)·M·S` and passes it in.

## Conventions & gotchas

- **`eleType` must match everywhere.** The mesh and every assembly call take
  `eleType`; a mismatch silently corrupts results. (`mesh.eleType` is also set
  by `generateColumnMesh`.)
- **`npg`** is Gauss points *per direction* (2 for Q4; use 3 for full Q8).
- **Always run from the project root** — paths in `main.m`/`runAllTests.m` are
  relative.
- **Figures:** `main.m` regenerates `tp3_tabla1/divergencia/fixedstress.png`.
  `tp3_esquema.png` is static geometry; its generator is kept commented at the
  end of `main.m`. The report pulls all four from `figs/` via `\graphicspath`.
- **The report is in Spanish.** Keep terminology/notation consistent with
  `report/Informe_TP3.tex` when editing docs or figure labels.
- **`old/` is legacy** single-field (displacement-only) elasticity. Reuse its
  B-matrix / Gauss-integration *style* if helpful, but it has no poroelastic
  coupling — don't wire it into the current pipeline.

## When changing things

- Editing `src/` → run `tests/runAllTests.m`.
- Changing mesh order, results, or `eleType` → re-run `main.m` to refresh
  `figs/`, then update the numbers/text in `report/Informe_TP3.tex`
  (iteration counts, node count, element-type wording).
- Derived scalars (`cv`, `p0`, `α²M/K_dr`, `ρ_fs`) depend only on material
  properties, not on the mesh — they don't change with Q4↔Q8.

# Static fragility audit — TP3 Terzaghi Poroelasticity

**Date**: 2026-05-30  
**Scope**: every file under `src/`, `meshes/`, plus `tp3.m`.  
**Status legend**:  ✅ guarded after E3 patch — ⚠ still open — 📐 design choice / known limitation

---

## 1. Hardcoded geometry assumptions

| Status | Location | Issue |
|---|---|---|
| ⚠ | `src/assembly/assembleLoad.m:24` | Top face detected as `y = max(mesh.nodes(:,2))`. Fails silently for any non-rectangular geometry. (Added `noTopFace` guard in E3 to catch the *empty* case, but does not detect *partial* top faces.) |
| ⚠ | `src/terzaghi_funs/buildCaseBC.m:9–14` | BC tests are literal `x = 0`, `x = W`, `y = 0`, `y = L`. Only valid for an axis-aligned box; non-box geometries silently get wrong BCs. |
| ⚠ | `src/postprocess/plotConsolidation.m:25–28` | `sort(yN)` over *all* nodes, then assumes the result describes a single column. Works only when `nx = 1`; with `nx > 1` produces visual noise (multiple FEM nodes per analytical-x). |

## 2. Dead code and duplication

| Status | Location | Issue |
|---|---|---|
| ⚠ | `src/terzaghi_funs/buildGlobalMatrices.m:5` | Constructs a *wrong* `C` matrix using a Lamé-by-Poisson formula, then immediately overwrites it on line 8 with the correct plane-strain matrix. Pure cruft. Slated for removal in Phase D2. |
| 📐 | `src/assembly/assembleMechanical.m` | `case 'Q4'` and `case 'Q8'` are byte-identical except for the string passed to `shapeFunctionsDer`. Kept per the explicit-template request (matches `old/getStiffnessMatrix.m`). |
| 📐 | `src/assembly/assemblePoroelastic.m` | Same observation. |

## 3. Silent failure paths

| Status | Location | Issue |
|---|---|---|
| ✅ | `src/elements/shapeFunctions.m`, `shapeFunctionsDer.m` | Added `otherwise → error('...:unknownEleType')` in E3. |
| ✅ | `src/assembly/assembleLoad.m` | Same. |
| ⚠ | `src/solvers/staggeredSolver.m:91–99` | On `maxIter` exhaustion the solver prints a warning but **continues with the non-converged state**, contaminating subsequent steps. No `status` flag is returned to the caller. |
| ⚠ | `src/solvers/fixedStressSolver.m` | Same: `iterHist(step) = maxIter` is set silently. |
| 📐 | `src/postprocess/analyticalTerzaghi.m:25` | Fixed 50-term truncation. For `t < L²/(c_v · 1000)` Gibbs phenomenon dominates. Limitation is documented in `developer_guide.md`; no run-time knob exposed. |

## 4. Floating-point tolerance assumptions

| Status | Location | Issue |
|---|---|---|
| ⚠ | `src/terzaghi_funs/buildCaseBC.m:3` | Hardcoded `tol = 1e-10`. For sub-millimetre `W`, mid-edge Q8 nodes at exact halves may *miss* the boundary. Recommend `tol = 1e-9 * max(W, L)`. |
| ⚠ | `src/assembly/assembleLoad.m:18` | Same `tol = 1e-10` — same scale-blindness. |

## 5. Mesh ↔ element-type coupling

| Status | Location | Issue |
|---|---|---|
| ⚠ | `meshes/generateColumnMesh.m` & `src/terzaghi_funs/buildGlobalMatrices.m` | `mesh.eleType` is not stored on the struct. Caller must pass a matching `eleType` to every assembly function — easy to mismatch after refactor (this exact mismatch caused the Q8 bug fixed earlier in the session). Slated for Phase E2. |
| ⚠ | `src/terzaghi_funs/buildGlobalMatrices.m` | Assembly assumes *all* elements share the same `eleType`. A mixed-type mesh would silently fail when the actual element has fewer/more nodes than expected. |

## 6. Numerical conditioning red flags

| Status | Location | Issue |
|---|---|---|
| ⚠ | `src/solvers/staggeredSolver.m:51` | Flow LU is built from `S/dt + H`. As `dt → 0`, `S/dt` dominates and the system becomes ill-conditioned. No warning issued; no condition-number check. |
| ⚠ | `src/solvers/fixedStressSolver.m:54` | Same observation, made worse by the added `Lstab/dt` term. |
| ⚠ | `tp3.m` & `fixedStressSolver` | The contraction factor ρ = α²M / (K_dr + α²M) is never printed; users cannot predict iteration count or set a realistic `tol`. |

## 7. Cosmetic / refactor

| Status | Location | Issue |
|---|---|---|
| 📐 | `tp3.m` | Two near-identical scenario blocks (Table 1 and Table 2). Slated for extraction into `runScenario` in Phase E5. |
| 📐 | `src/postprocess/analyticalTerzaghi.m` | Two near-identical loops compute `coeff_p` twice per time step. Vectorisation in Phase D3 will fold them into one matrix product. |

---

## Closure plan

Items marked ⚠ are tracked in the optimization plan:

| Plan phase | Closes |
|---|---|
| D1 / D2 | §2 dead code |
| D3 | §7 duplication in analytical solver |
| E1 | §6 (constitutive helper supports future modes) |
| E2 | §5 mesh ↔ eleType coupling |
| E5 | §7 scenario duplication |
| (not in current plan) | §1 hardcoded geometry, §3 solver `status` flag, §4 scale-aware tolerance, §6 conditioning warnings |

Items marked 📐 are deliberate design choices documented for users; no
patch is planned.

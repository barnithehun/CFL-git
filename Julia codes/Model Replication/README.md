# Model Replication

This folder holds the working copy of the structural model behind *Bankruptcy Resolution and the Macroeconomics of Cash Flow-Based Lending* (`Latex codes/Main/main5.tex`). The work has three steps:

1. Reproduce the paper's model results from a copy of `FinalVers/Baseline.jl`.
2. Cut computation time.
3. Recalibrate to the new targets from `Data Replication/calib_targets.R`.

The originals in `Julia codes/FinalVers/` and `Julia codes/Functions/` are read-only references. All model edits happen here.

## Rules for editing (read before changing any code)

The point of these rules is that the code keeps the structure it had in `FinalVers/`, so its author still recognises it after every change. Git records what changed, so edits need no extra markers or change log.

1. **Keep the structure by default.** Keep the same files, function names and signatures, and variable names (`SumPol` and its 18 columns, `phi_c`, `zeta_Rl`, `parameters()`, `gridsize()`, ...). Keep the same algorithm: Tauchen, then `DiscreteDP` with policy iteration, then the fixed-point loop on q, then `stat_dist`, then `sumSS`.
2. **Prefer changes inside a step.** A speed-up normally replaces lines inside an existing function, for example how `Q` is filled or how the q update loops. It does not move code into new modules, structs or files.
3. **Overhauls only when they are really worth it.** A larger restructuring is fine when the gain is large: a big speed-up, a correctness fix, or something the recalibration needs. Changes made only for elegance or flexibility are not worth it. That rules out renaming, reformatting, config structs in place of globals or `parameters()`, and generic wrappers.
4. **Ask before big changes.** Anything that rewrites more than a block, restructures files, or changes how the model is solved is proposed first, with its expected gain.
5. **Check results after each change.** `c_e`, `sumSS` and `sumSSsme` for the baseline and the reform must match the reference run, to rounding, unless the change is meant to move them. Record the reference values before the first speed change.

## Setup

The environment is `Project.toml` / `Manifest.toml` at the repo root, with Julia 1.12.7.

```julia
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))   # repo root (CFL-git)
Pkg.instantiate()
```

`Baseline.jl` is meant to be run block by block (for example with Shift+Enter in VS Code). Its first lines `cd` into this folder, so includes load the copies in `Functions/` and outputs are saved here.

**Memory:** `FirmOptim` builds a dense transition array `Q` of size n × m × n, with n = 1,059 states and m = 2,118 actions at the 46/23/46/46 grid. That is about 19 GB of RAM. `Abaseline.jl` avoids it and takes about 1.7 s per `FirmOptim` at the full grid.

## What is in this folder

| File | Origin | Changes |
|---|---|---|
| `Baseline.jl` | `FinalVers/Baseline.jl` at commit `6bd4edb` | Includes point to `Functions/` here, and a `cd(@__DIR__)` line was added |
| `Abaseline.jl` | Copy of `Baseline.jl` | The faster working version (option A). Only `FirmOptim` changes, in original lines 107–380, plus helper functions placed just before it. `FirmOptim` stores transition tables that depend on e only, in place of the dense `Q`. It solves the DP with `PFI_solve`, which follows the steps of QuantEcon's PFI, and runs the q update once per (e, a) pair. `PFI_solve` starts from the previous q iteration's values (warm start), so each solve needs fewer policy iterations. The q loop stops once k', b' and default are unchanged and q moved by at most 1e-10 (`q_tol`), instead of waiting for an exact repeat. `Q_sigma_mat` returns a sparse matrix, so the policy evaluation in `PFI_solve` uses a sparse solve. Everything after `FirmOptim` is unchanged |
| `Functions/*.jl`, `Functions/FCmodel/*.jl` | Every file that `Baseline.jl` includes | None |

On 2026-10-28, `Abaseline.jl` reproduced `Baseline.jl` bit for bit at k = b = 14, 20 and 30, for both the baseline and the reform. The comparison covered `SumPol`, `Fmat`, `c_e`, `sumSS` and `sumSSsme`. `Baseline.jl` cannot run at the full grid on a machine with about 16 GB of free RAM. There, `Abaseline.jl` reproduces the paper's untargeted small/large table and c_e = 654.3.

On 2026-10-29, the warm start was checked against the earlier `Abaseline.jl` results at the full grid. `SumPol`, `Fmat`, `c_e`, `sumSS` and `sumSSsme` were bit-identical for the baseline and the reform, and `FirmOptim` went from 5.8 s to 3.3 s per call. The tolerance stop cut the q loop from 15 to 10 iterations and the time to 2.3 s. Policies, `Fmat`, `sumSS` and `sumSSsme` stayed bit-identical; q, values and `c_e` moved by at most 1e-12. Sparse policy evaluation then cut the time to 1.7 s, with the same differences from the reference.

`Functions/PrintPolNew.jl` and `Functions/FCmodel/UpdateLoop.jl` were not copied because nothing uses them.

### Baseline.jl, block by block

Line numbers refer to the copy in this folder. This is the block order inherited from `FinalVers/`. It may later be rearranged to follow the paper's tables and figures, which counts as a big change under rule 4.

| Lines | Block | Paper output | Notes |
|---|---|---|---|
| 31–68 | `gridsize()`, `parameters()` | Calibrated parameter table | Values match main5 |
| 70–436 | `FirmOptim(wage, phi_c, zeta_Rl)` | — | The solver; main target of step 2 |
| 437–449 | Baseline (w = 1, θ = 0.3605, ζ = 2753.23) and reform (ζ = 1350, w = 1.0175) | — | The reform wage is hard-coded here |
| 452–459 | `sumSSsme`, `sumSS` | Targets, untargeted small/large and reform tables | |
| 462–479 | `Ushape`, `Xcross`, `QBplot`, `CapAlloc`, `DebtScedule` | Diagnostics, not in main5 | `DebtScedule` uses the stale `GridMake` |
| 481–500 | `dynsim2`, `plotPol`, `plotPDF`, `plotXE`, `plotTauDist` | Diagnostics | |
| 502–505 | `GamPol` | Appendix `GamPol.png` | |
| 507–511 | `improve.jl`, `plotPolHeatmaps` | Appendix `reform.png` | Including `improve.jl` also runs top-level code that saves `productivity_heatmap.png` and `combined_heatmaps.png` |
| 514–517 | `ZetaGam`, `ZetaGamPlot` | — | Broken: `SumPolZeta` is never created |
| 519–544 | GE wage by bisection (`FindWage`) | Reform wage | As saved, `zeta_Rl = 1350` is commented out and the bracket is (1.1, 1.3). Line 533 records w = 1.0175 for ζ/2 |
| 546–575 | Calibration: `ErrorFunc1` with Nelder–Mead, 8-hour limit | Calibrated parameters | Solves with `FirmOptim_Ext1`, a separate copy of the solver that uses Fcut 1200. `ErrorFunc2` is not defined |
| 579–602 | GE with shifted costs | Not in main5 | Uses old parameter values |
| 605–613 | `FindZeta` | — | Broken: called with 3 arguments, needs 4 |

`Functions/CapAlloc2.jl` is a script. It is run by hand after the baseline and reform blocks, and it makes the decomposition table from `SumPol`, `SumPol0`, `mu` and `mu0`.

### Functions

| File | Used for | Status |
|---|---|---|
| `FCmodel/fn_Tau_Q.jl` | Debt price q and CF share τ for each state–action pair | Core; calls `parameters()` twice per call |
| `FCmodel/EntryValue.jl` | Entry value `c_e` and entrant distribution `f0` | Core |
| `FCmodel/stat_dist.jl` | Stationary distribution `mu` | Core; reads the global `zeta_Rl` |
| `FCmodel/FindWage.jl` | Entry value at a given wage, for the GE bisection | Core |
| `FCmodel/ErrorFunc.jl` | Calibration objectives (`ErrorFunc1`, `_CF`, `_AB`) | Targets liquidation probability 0.55; main5 says 0.58 |
| `FCmodel/FirmOptim_Ext.jl` | Solver copies used by the calibration | Speed changes must be repeated here |
| `FCmodel/GridMake.jl` | Grids for `DebtScedule` and `ZetaGam` | Stale: unpacks `parameters()` in an old order |
| `sumSS.jl`, `sumSSsme.jl` | Aggregate and small/large moments | Each calls `stat_dist` again |
| `GamPol.jl` | `GamPol` (paper figure), `DebtScedule` (stale) | |
| `improve.jl` | `plotPolHeatmaps` | Runs code when included |
| `CapAlloc2.jl` | Decomposition table (script) | Hard-coded −0.003 / +0.003 adjustments (lines 107, 109) |
| `PolperSize.jl`, `CapAlloc.jl`, `dynsim2.jl`, `plotPol.jl`, `StatDist_plot.jl`, `PrintPolOld.jl` | Diagnostics and plots | |
| `dynsim.jl` | — | Included but never called |
| `ZetaGam.jl`, `ZetaGamPlot.jl` | — | Broken: 4-argument `FirmOptim` |
| `FindZeta.jl` | — | Ignores its argument and uses the global `zeta_Rs_alt` |

## Other FinalVers scripts (not copied yet)

| Script | Paper output | Calibration |
|---|---|---|
| `ABcase.jl` | AB-only column of the percentage-changes table | Older calibration (δ = 0.1, 32-point k/b grid) |
| `CFcase.jl` | CF-only column of the same table | Older calibration |
| `SmallCosts.jl` | SBRA table, targeted column | Same as Baseline |
| `Xallocation.jl` | Firm and capital allocation figures (4 runs) | Different calibration |
| `ShiftingCosts.jl`, `ShiftingFixedCosts.jl` | Not in main5 | Old calibration, broken |

No saved script reproduces the SBRA table's universal (0.93%) and control-reduction (θ 0.361 → 0.268) columns.

## Open questions about the original code

These are recorded here, not fixed. Each one is decided separately.

- **Firm law of motion.** `Fmat` moves P_x + P_D·γ of each state to exit, but P_D already includes P_x, and `stat_dist` also removes states that choose default. This may count defaults twice. Check it against the paper's equations.
- **Liquidation target.** `ErrorFunc1` uses 0.55, while main5's target table says 0.58.
- **Decomposition.** The ±0.003 in `CapAlloc2.jl` changes the policy margin from about 0.01 to 0.311 log points.
- **Size split.** main5 splits small and large firms by borrowing. `sumSSsme` uses k′ < Fcut, and `sumSS` uses k + x ≤ 2000.
- **Solver stop.** The q loop stops after 25 iterations with only a printed warning.

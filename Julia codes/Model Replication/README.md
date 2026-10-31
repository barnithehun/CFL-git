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

**Memory:** The original `FirmOptim` in `FinalVers/Baseline.jl` builds a dense transition array `Q` of size n × m × n, with n = 1,059 states and m = 2,118 actions at the 46/23/46/46 grid. That is about 19 GB of RAM. `Baseline.jl` here avoids it and takes about 0.65 s per `FirmOptim` at the full grid.

## What is in this folder

| File | Origin | Changes |
|---|---|---|
| `Baseline.jl` | `FinalVers/Baseline.jl` at commit `6bd4edb` | The faster working version. Includes point to `Functions/` here, and a `cd(@__DIR__)` line was added. Otherwise only `FirmOptim` changes, in original lines 107–380, plus helper functions placed just before it. `FirmOptim` stores transition tables that depend on e only, in place of the dense `Q`. It solves the DP with `PFI_solve`, which follows the steps of QuantEcon's PFI, and runs the q update once per (e, a) pair. `PFI_solve` starts from the previous q iteration's values (warm start), so each solve needs fewer policy iterations. The q loop stops once k', b' and default are unchanged and q moved by at most 1e-10 (`q_tol`), instead of waiting for an exact repeat. `Q_sigma_mat` returns a sparse matrix, so the policy evaluation in `PFI_solve` uses a sparse solve. Where x' falls on the x-grid for each (e', a) is computed once, before the q loop (`x_lo_ea`, `t_ea`, `interp_ea`), and the q update reads it from these tables. Everything after `FirmOptim` is unchanged |
| `baseline_calib.jl` | `Baseline.jl` up to the end of `FirmOptim` (the model), plus a calibration block and a model fit section | `parameters()` reads overrides (`par(:name, default)` on each line) and also returns `phi_c`, `zeta_Rl` and `DRS`, appended last so positional unpacking still works. Overriding `DRS` moves `alpha` and `nu` with it. `FirmOptim(wage, phi_c, zeta_Rl)` takes `phi_c` and `zeta_Rl` from `parameters()` when they are left out. It includes `FCmodel/Calibrate.jl` in place of `FCmodel/ErrorFunc.jl` and loads `CSV`. The baseline/reform results blocks of `Baseline.jl` are not copied; the model fit section replaces them |
| `Functions/*.jl`, `Functions/FCmodel/*.jl` | Every file that `Baseline.jl` includes | `FCmodel/fn_Tau_Q.jl` takes `beta` and `kappa` as arguments and builds no arrays; its callers in `FCmodel/FirmOptim_Ext.jl` and `GamPol.jl` pass them. `FCmodel/stat_dist.jl`, `sumSS.jl` and `sumSSsme.jl` take a keyword `zeta_R`, which defaults to the global `zeta_Rl` as before. `sumSS.jl` and `sumSSsme.jl` also take `rate_digits`, the digits the interest rate and spread are rounded to; the default 2 keeps their output unchanged, and `nothing` leaves them unrounded. `FCmodel/Calibrate.jl` is new. All other files unchanged |
| `calibration_outputs/` | Written by `calibrate` and `calibrate_thoroughly` | One `CALIB_<date>.csv` per calibration run; `calibrate_thoroughly` writes `CALIB_<date>_thorough.csv` and `CALIB_<date>_thorough_runs.csv` |

On 2026-10-28, the faster solver (then called `Abaseline.jl`) reproduced the original `Baseline.jl` bit for bit at k = b = 14, 20 and 30, for both the baseline and the reform. The comparison covered `SumPol`, `Fmat`, `c_e`, `sumSS` and `sumSSsme`. The original cannot run at the full grid on a machine with about 16 GB of free RAM. There, the faster solver reproduces the paper's untargeted small/large table and c_e = 654.3.

On 2026-10-29, the warm start was checked against the earlier results at the full grid. `SumPol`, `Fmat`, `c_e`, `sumSS` and `sumSSsme` were bit-identical for the baseline and the reform, and `FirmOptim` went from 5.8 s to 3.3 s per call. The tolerance stop cut the q loop from 15 to 10 iterations and the time to 2.3 s. Policies, `Fmat`, `sumSS` and `sumSSsme` stayed bit-identical; q, values and `c_e` moved by at most 1e-12. Sparse policy evaluation then cut the time to 1.7 s, with the same differences from the reference. Computing the x' grid positions once, and the new `fn_Tau_Q`, then cut it to 0.65 s with bit-identical results. The new `fn_Tau_Q` was also checked old against new in the original `Baseline.jl` and the three `FirmOptim_Ext` solvers at k = b = 14, and in `GamPol`: all bit-identical.

On 2026-10-30, the copy of the original `Baseline.jl` was dropped from this folder and `Abaseline.jl` was renamed `Baseline.jl`. The original is still in `FinalVers/Baseline.jl` and in git history.

On 2026-10-30, `baseline_calib.jl` reproduced the committed `Baseline.jl` (lines 1–570) bit for bit: `SumPol`, `Fmat`, `mu`, `c_e`, `sumSS` and `sumSSsme` for the baseline and the reform. `calib_fit` at the values in `parameters()` gives the same moments and c_e = 654.28. At k = b = 14, `FirmOptim_Ext1` and `FirmOptim` gave identical policies, `Fmat` and `sumSS`, and c_e within 1e-10, so the calibration can use `FirmOptim` itself.

On 2026-10-31, after the model fit section replaced the results blocks, `model_fit` at the values in `parameters()` reproduced the committed `Baseline.jl` bit for bit: the targeted `sumSS` moments, the small and large `sumSSsme` moments, the share of small firms and c_e. This used main5's targets with the share target set to the model's share at Fcut = 2000, and Fcut stayed at 2000. A recovery test at k = b = 20 (targets = the model's own moments at the current values, start 5–10% away) converged in 1,147 trials (390 s) to an error of 0.012. Every targeted moment came within 0.5% of its target, but the parameters came back 0–10% off (θ the most, P_x and c next). The objective is flat in places: policies are on a grid, and `sumSS` rounds the interest rate to two decimals. With unrounded rates the same test took 515 trials (178 s) and reached an error of 0.004, with every targeted moment within 0.26%. The parameters still came back 5–11% off (θ and ζ the most), so at this grid several parameter sets give nearly the same five moments.

The same test at the full grid, with unrounded targets, ran from two starts 5–10% from the truth in opposite directions. Each took about 660 trials and 11–12 minutes, two runs in parallel. Neither found the truth, where the error is 0. Nelder–Mead reported convergence at an error of 0.112 from start A and 0.053 from start B, with moments 1–3% off and parameters up to 19–21% off (P_x the most). The simplex had shrunk to points within about 1e-8 of each other in error. Moments jump where a policy switches grid points, so the error has local minima that are not the global one. The calibration code finds the right error at the truth; the weak point is Nelder–Mead on this rugged objective.

Restarting Nelder–Mead from where it stopped, with a fresh default simplex (each parameter +50%), was tested on both runs. From B, the first restart lowered the error from 0.053 to 0.038 (1,065 trials) and the second to 0.00006 (467 trials). The parameters were then within 0.01% of the truth, except `phi_c_hh` (−1%). A third restart gained nothing. From A, one restart (647 trials) gained only 3e-9: that run sits in a local minimum wider than the simplex. So restarts can finish a run that is close, but they cannot get a run out of a wrong basin.

`Functions/PrintPolNew.jl` and `Functions/FCmodel/UpdateLoop.jl` were not copied because nothing uses them.

### Calibration in baseline_calib.jl

The calibration block is the only place to edit. It sets four things:

- **`data_targets`**: a table with columns `moment`, `weighted`, `small` and `large`. With `targets_source = :R` it is read from `Data Replication/calib_targets.csv`, the last table of `calib_targets.R` plus its share of small firms w; run that script first to update it. With `:manual` it is `manual_targets`, main5's values.
- **`calib_moments`**: `(moment, row in sumSS, row in sumSSsme, error form)` for each moment. The error term is `|model − data| / data` (`:rel`), or `|model − data| / (1 − data)` (`:rel_1m`), which is `ErrorFunc1`'s AB share in place of the CF share. `:none` shows a moment without targeting it.
- **`calib_params`**: each internally calibrated parameter as `name = (start, lower, upper)`. Any number field of `parameters()` can be listed; a parameter that is not listed keeps its value from `parameters()`.
- **`valid`**: an optional joint condition on the trial values.

`calibrate` runs Nelder–Mead on the summed error terms. Each trial sets the trial values as overrides in `parameters()` through `with_params`, then runs `FirmOptim`, `EntryValue` and `sumSS`. The calibration and the model fit use unrounded interest rates and spreads (`rate_digits = nothing`). With the default two digits, the error stays flat until the rate crosses 0.005. A trial outside the bounds scores 1e10 instead of `Inf`: an `Inf` in the simplex makes Nelder–Mead's convergence measure NaN and stops it, and the default first simplex sets `phi_c_hh` to 1.23.

**Fcut is set after the fit.** It does not enter the solution while `pdef_exo_s = pdef_exo_l` (`model_fit` warns otherwise); it only splits small and large firms. The share of small firms (`sumSSsme`'s `share` row) changes only where Fcut passes a k′ value, so it is a step function that Nelder–Mead handles badly. `fit_Fcut` picks the split closest to the data's share of small firms. It keeps the current Fcut if that split is already the best, and otherwise uses the k′ value just above the split.

`calibrate` prints the lines to paste into `parameters()` and saves `calibration_outputs/CALIB_<date>.csv`. The file has one row per item, grouped by `section`: internally calibrated parameters (with their starts), Fcut, external parameters, grid sizes, targeted moments (data, model, error term), small and large firm moments, the share of small firms, and the total error, c_e and wage. The model fit section prints the same fit: total error, share of small firms, the targeted moments, and the small and large firm moments. `set_params!(CalRes.values)` makes the calibrated values (Fcut included) the active ones without pasting; `reset_params!()` undoes it.

**`calibrate_thoroughly`** runs only when `thorough = true` is set by hand in the calibration block; the default is `calibrate`. It adds restarts and several starts, for the two failures found in the recovery tests above:

- From each start it runs Nelder–Mead, then restarts it from where it stopped, until a restart improves the error by less than `restart_tol` (1e-6), or after `max_restarts` (5).
- The default starts are the start in `calib_params` and two variants that move each start by ±`variant_pct` (10%) with alternating signs (+,−,+,… and −,+,−,…). A start outside the bounds or failing `valid` is skipped with a warning.
- `starts = [CalRes.values]` refines an earlier result instead; fields that are not in `calib_params`, such as Fcut, are ignored.
- It keeps the best run. The first run from the start in `calib_params` is the same as `calibrate`'s, so the result is never worse.

It saves the same CSV as `calibrate`, named `CALIB_<date>_thorough.csv`, with the start of the best run in the start column. `CALIB_<date>_thorough_runs.csv` lists every run: start, restart, error, gain, trials, minutes, the values and the best run. If the starts end at different values with similar errors, the moments do not pin down the parameters. It takes several times as long as `calibrate`; `calibrate` itself is unchanged. Tested on 2026-10-31 at k = b = 20 with 15 iterations per run: the first run gave exactly `calibrate`'s error, the best run was the one kept, restarts stopped at a zero gain, `starts = [CalRes.values]` worked, and a start outside the bounds was skipped.

### Baseline.jl, block by block

Line numbers refer to the copy in this folder. This is the block order inherited from `FinalVers/`. It may later be rearranged to follow the paper's tables and figures, which counts as a big change under rule 4.

| Lines | Block | Paper output | Notes |
|---|---|---|---|
| 31–68 | `gridsize()`, `parameters()` | Calibrated parameter table | Values match main5 |
| 70–192 | `Qv!`, `Q_sigma_mat`, `bellman!`, `PFI_solve` | — | Helpers for `FirmOptim`; not in the original |
| 194–544 | `FirmOptim(wage, phi_c, zeta_Rl)` | — | The solver; main target of step 2 |
| 546–558 | Baseline (w = 1, θ = 0.3605, ζ = 2753.23) and reform (ζ = 1350, w = 1.0175) | — | The reform wage is hard-coded here |
| 561–568 | `sumSSsme`, `sumSS` | Targets, untargeted small/large and reform tables | |
| 571–588 | `Ushape`, `Xcross`, `QBplot`, `CapAlloc`, `DebtScedule` | Diagnostics, not in main5 | `DebtScedule` uses the stale `GridMake` |
| 590–609 | `dynsim2`, `plotPol`, `plotPDF`, `plotXE`, `plotTauDist` | Diagnostics | |
| 611–614 | `GamPol` | Appendix `GamPol.png` | |
| 616–620 | `improve.jl`, `plotPolHeatmaps` | Appendix `reform.png` | Including `improve.jl` also runs top-level code that saves `productivity_heatmap.png` and `combined_heatmaps.png` |
| 623–626 | `ZetaGam`, `ZetaGamPlot` | — | Broken: `SumPolZeta` is never created |
| 628–653 | GE wage by bisection (`FindWage`) | Reform wage | As saved, `zeta_Rl = 1350` is commented out and the bracket is (1.1, 1.3). Line 642 records w = 1.0175 for ζ/2 |
| 655–684 | Calibration: `ErrorFunc1` with Nelder–Mead, 8-hour limit | Calibrated parameters | Solves with `FirmOptim_Ext1`, a separate copy of the solver that uses Fcut 1200. `ErrorFunc2` is not defined |
| 688–711 | GE with shifted costs | Not in main5 | Uses old parameter values |
| 714–721 | `FindZeta` | — | Broken: called with 3 arguments, needs 4 |

`Functions/CapAlloc2.jl` is a script. It is run by hand after the baseline and reform blocks, and it makes the decomposition table from `SumPol`, `SumPol0`, `mu` and `mu0`.

### Functions

| File | Used for | Status |
|---|---|---|
| `FCmodel/fn_Tau_Q.jl` | Debt price q and CF share τ for each state–action pair | Core; takes `beta` and `kappa` from the caller |
| `FCmodel/EntryValue.jl` | Entry value `c_e` and entrant distribution `f0` | Core |
| `FCmodel/stat_dist.jl` | Stationary distribution `mu` | Core; reads the global `zeta_Rl` unless `zeta_R` is passed |
| `FCmodel/Calibrate.jl` | `calibrate`, `calibrate_thoroughly`, `default_starts`, `model_fit`, `print_fit`, `fit_Fcut`, `read_targets`, `save_calibration`, `calib_fit`, `calib_error`, and the overrides `set_params!`, `reset_params!`, `with_params` | Used by `baseline_calib.jl` only |
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
- **Liquidation target.** `ErrorFunc1` uses 0.55, while main5's target table says 0.58. `baseline_calib.jl` uses 0.58.
- **Decomposition.** The ±0.003 in `CapAlloc2.jl` changes the policy margin from about 0.01 to 0.311 log points.
- **Size split.** main5 splits small and large firms by borrowing. `sumSSsme` uses k′ < Fcut, and `sumSS` uses k + x ≤ 2000.
- **Solver stop.** The q loop stops after 25 iterations with only a printed warning.

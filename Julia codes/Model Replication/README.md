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
| `baseline_calib.jl` | `Baseline.jl` up to the end of `FirmOptim` (the model), plus a calibration block and a model fit section | `parameters()` reads overrides (`par(:name, default)` on each line) and also returns `phi_c`, `zeta_Rl` and `DRS`, appended last so positional unpacking still works. Overriding `DRS` moves `alpha` and `nu` with it. `FirmOptim(wage, phi_c, zeta_Rl)` takes `phi_c` and `zeta_Rl` from `parameters()` when they are left out. It includes `FCmodel/Calibrate.jl` in place of `FCmodel/ErrorFunc.jl` and loads `CSV`. Three switches for the results experiments are appended to `parameters()`: `zeta_Rs` and `phi_c_s`, ζ and θ for small firms (k′ < Fcut), which follow `zeta_Rl` and `phi_c` unless set, and `perf_cred` (default 0). In `FirmOptim`, small firms' ζ and θ follow its `zeta_Rl` and `phi_c` arguments unless they are set as overrides, and `perf_cred = 1` allows negative dividends. With `perf_cred = 1`, `PFI_solve` also starts cold in each q iteration, as the original solver did: some debt choices then tie, and the warm start picked a different tied choice in each q iteration, so the q loop cycled. At their defaults the model is bit-identical to before. The baseline/reform results blocks of `Baseline.jl` are not copied; the RESULTS section after the model fit replaces them (see below) |
| `MiscResults.jl` | The diagnostics and experiments of `Baseline.jl`, `SmallCosts.jl`, `ShiftingCosts.jl` and `improve.jl` that the paper does not use | New. Labelled blocks M1–M11; see below |
| `Functions/*.jl`, `Functions/FCmodel/*.jl` | Every file that `Baseline.jl` includes | `FCmodel/fn_Tau_Q.jl` takes `beta` and `kappa` as arguments and builds no arrays; its callers in `FCmodel/FirmOptim_Ext.jl` and `GamPol.jl` pass them. `FCmodel/stat_dist.jl`, `sumSS.jl` and `sumSSsme.jl` take a keyword `zeta_R`, which defaults to the global `zeta_Rl` as before. `sumSS.jl` and `sumSSsme.jl` also take `rate_digits`, the digits the interest rate and spread are rounded to; the default 2 keeps their output unchanged, and `nothing` leaves them unrounded. Their debt to collateral ratio is b′/k′, as in the data (debt over PPE, inventories and receivables), where the original divided by k′ + x; ratios above 3 are left out, as `calib_targets.R` does (`d2c_cut`). `CapAlloc2.jl` is now a function, `CapAlloc2(B, R)`, with the script's body. `FCmodel/Calibrate.jl`, `Results.jl` and `Xallocation.jl` (the results block of `FinalVers/Xallocation.jl` as a function) are new. All other files unchanged |
| `calibration_outputs/` | Written by `calibrate` and `calibrate_thoroughly` | One `CALIB_<date>.csv` per calibration run; `calibrate_thoroughly` writes `CALIB_<date>_thorough.csv`, `CALIB_<date>_thorough_runs.csv` and `CALIB_<date>_thorough_screen.csv`. The AB-only and CF-only recalibrations write `CALIB_<date>_AB.csv` and `CALIB_<date>_CF.csv` |
| `results_outputs/` | Written by the RESULTS section and `MiscResults.jl` | `<date>/RESULTS.csv` and the paper figures per results run; `misc_<date>/` holds `MiscResults.jl`'s figures |

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

- **`data_targets`**: a table with columns `moment`, `weighted`, `small` and `large`. It is read from `Data Replication/calib_targets.csv`, the last table of `calib_targets.R` plus its share of small firms w; run that script first to update it.
- **`calib_moments`**: `(moment, row in sumSS, row in sumSSsme, targeted)` for each moment. A targeted moment's error term is `|model − data| / data`; `targeted = false` shows a moment without targeting it.
- **`calib_params`**: each internally calibrated parameter as `name = (start, lower, upper)`. Any number field of `parameters()` can be listed; a parameter that is not listed keeps its value from `parameters()`.
- **`valid`**: an optional joint condition on the trial values.

`calibrate` runs Nelder–Mead on the summed error terms. Each trial sets the trial values as overrides in `parameters()` through `with_params`, then runs `FirmOptim`, `EntryValue` and `sumSS`. The calibration and the model fit use unrounded interest rates and spreads (`rate_digits = nothing`). With the default two digits, the error stays flat until the rate crosses 0.005. A trial outside the bounds scores 1e10 instead of `Inf`: an `Inf` in the simplex makes Nelder–Mead's convergence measure NaN and stops it, and the default first simplex sets `phi_c_hh` to 1.23.

**Fcut is set after the fit.** It does not enter the solution while `pdef_exo_s = pdef_exo_l` (`model_fit` warns otherwise); it only splits small and large firms. The share of small firms (`sumSSsme`'s `share` row) changes only where Fcut passes a k′ value, so it is a step function that Nelder–Mead handles badly. `fit_Fcut` picks the split closest to the data's share of small firms. It keeps the current Fcut if that split is already the best, and otherwise uses the k′ value just above the split.

`calibrate` prints the lines to paste into `parameters()` and saves `calibration_outputs/CALIB_<date>.csv`. The file has one row per item, grouped by `section`: internally calibrated parameters (with their starts), Fcut, external parameters, grid sizes, targeted moments (data, model, error term), small and large firm moments, the share of small firms, and the total error, c_e and wage. The model fit section prints the same fit: total error, share of small firms, the targeted moments, and the small and large firm moments. `set_params!(CalRes.values)` makes the calibrated values (Fcut included) the active ones without pasting; `reset_params!()` undoes it.

**`calibrate_thoroughly`** runs only when `thorough = true` is set by hand in the calibration block; the default is `calibrate`, which is unchanged. On 2026-10-31, at the full grid with the R targets, `calibrate` stopped at an error of 0.188. The first version of `calibrate_thoroughly` (restarts from the start and two ±10% or ±20% variants) reached 0.096 and 0.066, in 87 and 117 minutes. Nearly all the gain came from the variant starts, not the restarts, and each best result came from a different start. So the current version spends its time on starts:

- **Screening.** It computes the error at `n_screen` (100) points spread evenly over `screen_box` (a Halton sequence), plus the start in `calib_params`, one model solve each. The best `n_starts` (4) become the starts. `screen_box` gives `name = (low, high)`; a parameter left out gets start × (0.5, 1.5), cut to its bounds. Points outside the bounds or failing `valid` are dropped. `starts = [CalRes.values]` skips the screening and refines an earlier result; fields not in `calib_params`, such as Fcut, are ignored.
- **Restarts.** From each start it runs Nelder–Mead, then restarts it from where it stopped, until a restart improves the error by less than `restart_tol` (1e-6), or after `max_restarts` (3). It keeps the best run over all starts.
- **Transformed parameters** (`transform = true`). Nelder–Mead searches the logit of a parameter with two finite bounds and log(p − lower) of one with only a lower bound, so no trial leaves the bounds. The first simplex moves each of these by `simplex_step` (0.4: ×1.5 for a log, 0.80 → 0.86 for a logit). With `transform = false` it searches the parameters themselves, as `calibrate` does.
- **Workers.** The screening and the starts run on `n_workers` (= `n_starts`) separate Julia processes. Threads would not work, because the overrides in `parameters()` are one global Dict. Each worker loads the model part of `baseline_calib.jl` (everything above `########## CALIBRATION`), so save the file before running. Overrides set with `set_params!` before the call are passed on. `n_workers = 0` runs everything in the current process.

It saves the same CSV as `calibrate`, named `CALIB_<date>_thorough.csv`, with the start of the best run in the start column. `CALIB_<date>_thorough_runs.csv` lists every run: start, restart, error, gain, trials, minutes, the values and the best run. `CALIB_<date>_thorough_screen.csv` lists every screened point, sorted by error. If the starts end at different values with similar errors, the moments do not pin down the parameters. Tested on 2026-11-01 at k = b = 20 with 10 iterations per run and 2 workers:
- the screening ranked 9 points, and a worker's error at the start matched the main process's;
- the best of 2 starts was kept, and all three CSVs were written;
- `starts = [CalRes.values]` with `n_workers = 0` and `transform = false` worked.

The whole first run took 1.5 minutes.

### Results in baseline_calib.jl

The RESULTS section after the model fit makes the paper's model results at the calibrated values: `CalRes.values` from the calibration in the session, or a saved calibration through `CalVals = load_calibration("calibration_outputs/CALIB_<date>.csv")`. `CalVals = NamedTuple()` uses the values written in `parameters()`. The helpers are in `Functions/Results.jl`.

Each equilibrium is a set of overrides in `parameters()` plus a wage. `solve_eq(vals; wage)` solves one and returns its overrides and wage, `SumPol`, `Fmat`, `f0`, `mu`, `c_e`, and `sumSS` and `sumSSsme` with unrounded rates. `ge_eq(vals, c_e)` finds the wage at which free entry holds with the baseline entry cost: bisection on the wage with `FindWage`, as in `Baseline.jl`'s GE block, until the bracket is narrower than 1e-6. The default bracket is (0.9, 1.3), and it stops with an error if the bracket does not contain the wage. This replaces the reform's hard-coded wage of 1.0175. The entry value can jump with the wage, because policies are on a grid and q and the liquidation decision jump with them. Then no wage gives free entry exactly, and the bisection ends at the jump. `ge_eq` keeps the end of the last bracket that is closer to free entry and warns when the gap is above 0.1% of c_e.

| Block | Equilibria | Paper output |
|---|---|---|
| Baseline and reform | `BaseEq` (w = 1); `ReformEq`: ζ/2, GE wage | tab:SSvalues_firms (`size_table`), tab:SSvalues_GE (`ge_table`); tab:targets and tab:untargeted are `Fit` |
| Decomposition | the same | tab:decomposition (`CapAlloc2`, without the ±0.003) |
| Allocation | `Reform10Eq`: ζ/10; `PerfEq`: `perf_cred = 1`, bracket (1, 1.6); both GE | fig:firmalloc, fig:capalloc (`Xallocation`); the perfect-credit change in Y/L is main5's footnote |
| GamPol | `BaseEq` | appendix GamPol.png |
| SBRA | `TargetedEq`: `zeta_Rs` from `find_zeta_Rs`; `UniversalEq`: the same ζ for all firms; `ControlEq`: targeted, and θ × ζ_s/ζ for all firms; all GE | tab:percentage_changes2 (`change_table`) |
| AB only, CF only | θ = 0 (`phi_c`) or φ_a = 0 (`phi_af`); recalibrated, then ζ/2 with GE wage | tab:percentage_changes, columns 2 and 3 |

- **`find_zeta_Rs`** searches a grid of 40–100% of ζ, at the baseline wage as in `SmallCosts.jl`, for the small-firm cost that lowers small firms' liquidation probability by 13 pp (main5 section 6). `SmallCosts.jl` noted that bisection fails here because the probability is not monotone in ζ_s.
- **The control reduction** cuts θ in proportion to the targeted cut, as main5 describes (0.361 → 0.268). The cut applies to all firms, as in the old code; `phi_c_s` would limit it to small firms.
- **AB only and CF only:** `calibrate` recalibrates pc, P_x and ζ to the debt to collateral ratio, the interest rate and the liquidation probability (main5), with the case's fixed value as an outer override (`with_params`), from the calibrated values. It needs `data_targets`, `calib_moments` and `options` from the calibration block, and saves `CALIB_<date>_AB.csv` and `_CF.csv`. With `recalibrate_ABCF = false`, `AB_file` and `CF_file` point to saved files instead.
- **Saving:** `save_results` writes `results_outputs/<date>/RESULTS.csv`, with one row per table cell and per equilibrium (its overrides, wage and c_e), and the figures as PNGs. Rates are unrounded in the file and rounded to 3 digits when printed.

Running `baseline_calib.jl` as a whole (for example with `include`) now also runs the RESULTS section after the calibration, including the AB/CF recalibrations. A GE solve takes about 15 s at the full grid. Solves in which the q loop hits its 25-iteration cap (some steps of perfect credit and of the control reduction) take several times longer.

### MiscResults.jl

`MiscResults.jl` holds the results that the paper does not use, in labelled blocks that name their source:

| Block | Content | From |
|---|---|---|
| M1 | `Ushape`, `Xcross`, `QBplot`, `CapAlloc` | `Baseline.jl` lines 571–580 |
| M2 | `DebtScedule`: commented out, it follows an older model | `Baseline.jl` lines 583–588 |
| M3 | `dynsim2`, `plotPol` | `Baseline.jl` lines 590–599 |
| M4 | `plotPDF`, `plotCDF`, `plotXE`, `plotTauDist`, `ProdDist` | `Baseline.jl` lines 601–609 |
| M5 | `plotPolHeatmaps` (appendix reform.png) | `Baseline.jl` lines 616–620 |
| M6 | `CapAlloc2`'s plots; the decomposition with the ±0.003 | `CapAlloc2.jl` |
| M7 | Policies over ζ = 0…8000 (`ZetaGamPlot`), solved with `FirmOptim` in place of the broken `ZetaGam` | `Baseline.jl` lines 623–626 |
| M8 | Targeted vs universal cuts at the baseline wage, ζ_s = 10–100% of ζ (smallcosts.png) | `SmallCosts.jl` lines 545–653 |
| M9 | Heatmaps over small and large firms' ζ, at the baseline wage | `improve.jl` lines 1–84 (commented out there) |
| M10 | ζ/2 and θ/2, for all firms and for small firms only, GE | `Baseline.jl` lines 688–711, `ShiftingCosts.jl` |
| M11 | θ − shift at ζ/2, at the baseline wage | `improve.jl` lines 153–210 (commented out there) |

It runs three ways. `run_misc = true` at the end of the RESULTS section includes it once, on `BaseEq` and `ReformEq`. On its own, run `baseline_calib.jl` up to the CALIBRATION block (the model) and include `MiscResults.jl`: its SETUP block solves the baseline and the reform again from the newest `results_outputs/<date>/RESULTS.csv`, or from `results_file` if set, in a few seconds. After SETUP, each block runs on its own. SETUP sets the globals the original functions read (`SumPol`, `Fmat`, `mu`, `SumPol0`, `mu0`, ...). Figures go to `results_outputs/misc_<date>/`.

### Checks of the results code (2026-10-31)

All at the values in `parameters()` (main5's), with outputs written outside the repo.

- **The three switches** (`zeta_Rs`, `phi_c_s`, `perf_cred`), at their defaults and also set explicitly to them: `SumPol`, `Fmat`, `c_e`, `sumSS` and `sumSSsme` are bit-identical to `baseline_calib.jl` before the change, at the full grid, for the baseline and for the reform at w = 1.0175. Rechecked after the perfect-credit cold start.
- **Against the FinalVers originals**, at k = b = 14, since the originals build the dense Q. Each original is `FinalVers/Baseline.jl`'s model with the named script's `FirmOptim`.
  - `Baseline.jl`, baseline and reform at w = 1.0175: policies identical; `Fmat`, `mu`, `sumSS` and `sumSSsme` bit-identical; `SumPol` within 4e-7 and `c_e` within 2e-9 (the tolerance stop).
  - `CapAlloc2.jl`: `CapAlloc2(...; adjust = 0.003)` is bit-identical to the script.
  - GE wage for ζ/2: `ge_eq` gives 1.0086929, `Baseline.jl`'s bisection on `FindWage` gives 1.0086926 (within the 1e-6 stop). That took 688 s, `ge_eq` 3.5 s.
  - `SmallCosts.jl`'s `FirmOptim(w, θ, ζ, ζ_s)` against the `zeta_Rs` override, targeted (ζ_s = 2030) and universal: the same match as the baseline.
  - `ABcase.jl` (θ = 0) and `CFcase.jl` (φ_a = 0), their `FirmOptim` with the older interpolation code: the same match as the baseline.
  - `Xallocation.jl`'s perfect credit: at w = 1.2 the same match as the baseline. At w = 1.1, `c_e` agrees to 1e-12 but some debt choices differ, because tied debt choices are broken differently. Under perfect credit, values, c_e and the wage are pinned down, and debt moments depend on how ties are broken (up to 0.005 here).
- **The whole RESULTS section and `MiscResults.jl`** ran end to end at k = b = 20, with the AB/CF recalibrations cut to 3 iterations, and wrote every table and figure. `MiscResults.jl` run on its own in a fresh session solved the baseline and the reform again from `RESULTS.csv`, with the same c_e.
- **Against main5's tables, at the full grid and main5's inputs** (`solve_eq` at w = 1.0175 for ζ/2, and SmallCosts.jl's ζ_s = 2030 and w = 1.0081):
  - tab:decomposition with `adjust = 0.003` is main5's to the digit: 1.735, 0.295, 2.802, −1.672, 0.311.
  - The reform column of tab:percentage_changes matches to main5's digits (mass 7.10%, Y/L 1.75%, liquidation −17.71 pp, debt/collateral 2.57 pp, rate −0.34 pp, CF reliance 14.7 pp, CF share 9.19 pp), except K/L: 2.65% here, while main5 reports 1.34% in this table and 2.92% in tab:SSvalues_GE.
  - tab:SSvalues_GE matches (mass 10.715, 933 employees, K/L 3.52, Y/L +1.75%).
  - The reform column of tab:SSvalues_firms is close: small firms' liquidation 0.634 against main5's 0.63, CF share 0.78 against 0.72.
  - The SBRA targeted column matches only in Y/L and the wage (0.81%). Liquidation is −5.4 pp against main5's −8.23, and CF reliance −5.1 pp against +5.3. The k = 14 check above shows the new code matches SmallCosts.jl's solver, so main5's column probably comes from an older run.
- **The full RESULTS section and `MiscResults.jl` at the full grid** ran end to end at main5's values (the AB/CF recalibrations cut to 3 iterations). `ge_eq` gave: ζ/2 at w = 1.01755, the right side of a jump (see below), so mass +6.6% instead of +7.1%; ζ/10 at 1.0260; perfect credit at 1.1455, with Y/L +14.6% (main5's footnote: about 14%); targeted ζ_s = 2037 at 1.0081 (Y/L +0.81%); universal at 1.0090 (Y/L +0.90%, main5 0.93%); control reduction (θ = 0.267) at 0.9924 (Y/L −0.76%, main5 −0.54%).
- **The entry value jumps with the wage.** At k = b = 20, in the three SBRA experiments, entry value − c_e jumps from about +8 to about −19 between wages 1e-6 apart, so no wage gives free entry exactly; the universal case shows it too, and it uses none of the switches. At the full grid, the ζ/2 reform jumps from +9.5 at w = 1.0175491 to −6.3 at 1.0175499; main5's 1.0175 is on the left side (c_e 663.83 against 654.28), and `ge_eq` keeps the right side, the one closer to free entry. The CF-only reform also jumps (+19.7 to −7.4 at w = 1.01221). Decided on 2026-11-01: the results report the side closer to free entry, as `ge_eq` does.

### Baseline.jl, block by block

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

`Functions/CapAlloc2.jl` was a script run by hand after the baseline and reform blocks, on the globals `SumPol`, `SumPol0`, `mu` and `mu0`. It is now the function `CapAlloc2(B, R)`. To use it after `Baseline.jl`'s blocks: `CapAlloc2((gridsize = gridsize(), params = parameters(), SumPol = SumPol, mu = mu), (SumPol = SumPol0, mu = mu0); adjust = 0.003)` gives the script's numbers.

### Functions

| File | Used for | Status |
|---|---|---|
| `FCmodel/fn_Tau_Q.jl` | Debt price q and CF share τ for each state–action pair | Core; takes `beta` and `kappa` from the caller |
| `FCmodel/EntryValue.jl` | Entry value `c_e` and entrant distribution `f0` | Core |
| `FCmodel/stat_dist.jl` | Stationary distribution `mu` | Core; reads the global `zeta_Rl` unless `zeta_R` is passed |
| `FCmodel/Calibrate.jl` | `calibrate`, `calibrate_thoroughly` (with `thorough_start`, `thorough_error`, `start_workers`, `halton`, `to_z`/`from_z`, `StepSimplexer`, `default_box`), `model_fit`, `print_fit`, `fit_Fcut`, `read_targets`, `save_calibration`, `calib_fit`, `calib_error`, and the overrides `set_params!`, `reset_params!`, `with_params` | Used by `baseline_calib.jl` only |
| `FCmodel/FindWage.jl` | Entry value at a given wage, for the GE bisection | Core |
| `FCmodel/ErrorFunc.jl` | Calibration objectives (`ErrorFunc1`, `_CF`, `_AB`) | Targets liquidation probability 0.55; main5 says 0.58 |
| `FCmodel/FirmOptim_Ext.jl` | Solver copies used by the calibration | Speed changes must be repeated here |
| `FCmodel/GridMake.jl` | Grids for `DebtScedule` and `ZetaGam` | Stale: unpacks `parameters()` in an old order |
| `Results.jl` | `solve_eq`, `ge_eq`, `find_zeta_Rs`, the tables (`size_table`, `ge_table`, `change_table`, `print_table`), `save_results`, `load_eqs`, `newest_results`, `load_calibration` | Used by the RESULTS section and `MiscResults.jl` |
| `sumSS.jl`, `sumSSsme.jl` | Aggregate and small/large moments | Each calls `stat_dist` again |
| `GamPol.jl` | `GamPol` (paper figure), `DebtScedule` (stale) | `DebtScedule` follows an older model version (the lender pays ζ in `Pi_reo`, old liquidation rule) |
| `improve.jl` | `plotPolHeatmaps` | Including it runs `plotPolHeatmaps(SumPol, SumPol0, mu, 1)` on the globals; its other code is commented out |
| `CapAlloc2.jl` | `CapAlloc2(B, R; adjust = 0.0)`: decomposition table and the script's plots | `adjust = 0.003` gives the script's −0.003 / +0.003 adjustments |
| `Xallocation.jl` | `Xallocation(runs, perf)`: firm and capital allocation figures | The results block of `FinalVers/Xallocation.jl`, on equilibria solved elsewhere |
| `PolperSize.jl`, `CapAlloc.jl`, `dynsim2.jl`, `plotPol.jl`, `StatDist_plot.jl`, `PrintPolOld.jl` | Diagnostics and plots | |
| `dynsim.jl` | — | Included but never called |
| `ZetaGam.jl`, `ZetaGamPlot.jl` | — | Broken: 4-argument `FirmOptim` |
| `FindZeta.jl` | — | Ignores its argument and uses the global `zeta_Rs_alt` |

## Other FinalVers scripts

None of these is copied; what they compute is now in `baseline_calib.jl`'s RESULTS section or in `MiscResults.jl`.

| Script | Paper output | Calibration | Now |
|---|---|---|---|
| `ABcase.jl` | AB-only column of the percentage-changes table | Older calibration (δ = 0.1, 32-point k/b grid) | RESULTS: θ = 0, recalibrated |
| `CFcase.jl` | CF-only column of the same table | Older calibration. Its calibration (`ErrorFunc_CF`, `FirmOptim_Ext_CF`) solved with φ_a = 0.4 and `nul_e` = 10 hard-coded, so it did not calibrate the CF-only model | RESULTS: φ_a = 0 (`phi_af`), recalibrated |
| `SmallCosts.jl` | SBRA table, targeted column | Same as Baseline | RESULTS (`zeta_Rs`); its targeted vs universal sweep is `MiscResults.jl` M8 |
| `Xallocation.jl` | Firm and capital allocation figures (4 runs) | Different calibration | RESULTS (`perf_cred`, `Xallocation`) |
| `ShiftingCosts.jl` | Not in main5 | Old calibration | `MiscResults.jl` M10 (`phi_c_s`) |
| `ShiftingFixedCosts.jl` | Not in main5 | Old calibration | Not ported: the lender pays a fixed cost in `Pi_reo`, a model change |

No saved script reproduced the SBRA table's universal (0.93%) and control-reduction (θ 0.361 → 0.268) columns; the RESULTS section computes both.

## Open questions about the original code

These are recorded here, not fixed. Each one is decided separately.

- **Firm law of motion.** `Fmat` moves P_x + P_D·γ of each state to exit, but P_D already includes P_x, and `stat_dist` also removes states that choose default. This may count defaults twice. Check it against the paper's equations.
- **Liquidation target.** `ErrorFunc1` uses 0.55, while main5's target table says 0.58. `baseline_calib.jl` takes its targets from `calib_targets.R` instead.
- **Decomposition.** The ±0.003 in `CapAlloc2.jl` changes the policy margin from about 0.01 to 0.311 log points. The RESULTS section reports the margins without it (`adjust = 0`); `MiscResults.jl` M6 prints both.
- **Size split.** main5 splits small and large firms by borrowing. `sumSSsme` uses k′ < Fcut, and `sumSS` uses k + x ≤ 2000.
- **Solver stop.** The q loop stops after 25 iterations with only a printed warning.

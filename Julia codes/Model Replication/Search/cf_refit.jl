# cf_refit.jl <run id> [name=value ...]: a thorough refit of the CF-only economy (no asset-based lending, phi_af = 0) of a run,
# and its zeta/2 reform. Search.jl fits the CF-only economy with calibrate_lazy (f_reltol 0.002) from the full model's values;
# here calibrate_fully (100-point screen, 4 starts, tighter tolerance) fits the same three parameters (pc, pdef_exo_l, zeta_Rl)
# to the same three aggregate targets (debt to collateral, interest rate, liquidation probability). Also re-solves the reform
# at the run's lazy CF-only values (runs/<id>/CALIB_CF.csv) as a check. Writes runs/<id>/CALIB_CF_full.csv and
# better_off/<id>_<grid>_cfrefit.csv.
using Dates, CSV, DataFrames, Printf
const SEARCH_DIR = @__DIR__
const MODEL_DIR  = normpath(joinpath(SEARCH_DIR, ".."))
const MODEL_FILE = joinpath(MODEL_DIR, "Baseline_noise.jl")
const DATA_DIR   = normpath(joinpath(MODEL_DIR, "..", "..", "Data Replication"))
ENV["GKSwstype"] = "100"
using LinearAlgebra
BLAS.set_num_threads(1)
let code = read(MODEL_FILE, String)
    i = findfirst("########## " * "CALIBRATION", code)
    task_local_storage(() -> include_string(Main, code[1:prevind(code, first(i))], MODEL_FILE), :SOURCE_PATH, MODEL_FILE)
end
cd(MODEL_DIR)
include(joinpath(MODEL_DIR, "Functions", "Results.jl"))

id = ARGS[1]
grid_over = Dict{Symbol, Float64}(Symbol(split(a, "=")[1]) => parse(Float64, split(a, "=")[2]) for a in ARGS[2:end])
run_dir = joinpath(SEARCH_DIR, "runs", id)
diag = first(CSV.read(joinpath(run_dir, "DIAG.csv"), DataFrame; types = String))
dv(k) = (k in names(diag) && !ismissing(diag[k]) && diag[k] != "") ? parse(Float64, diag[k]) : nothing

# external settings as in Search.jl: the recorded values, then the run's set: overrides (DIAG columns x_<name>)
ext = Dict{Symbol, Float64}()
for k in ("beta", "delta", "phi_a", "phi_af", "kappa", "rho_e", "sigma_e", "nul_e", "DRS", "size_by_debt", "d2c_cut", "tauchen_lo")
    v = dv(k); v === nothing || (ext[Symbol(k)] = v)
end
for k in names(diag)
    startswith(k, "x_") && dv(k) !== nothing && (ext[Symbol(k[3:end])] = dv(k))
end
merge!(ext, grid_over)
vals = load_calibration(joinpath(run_dir, "CALIB.csv"))
reset_params!(); set_params!((; ext...)); set_params!(vals)
g = gridsize()
gtag = "$(g.x_size)_$(g.e_size)_$(g.k_size)_$(g.b_size)"
println("$id on grid $gtag; externals ", ext); flush(stdout)

tag = diag["targets"]
data = read_targets(joinpath(DATA_DIR, "calib_targets_$tag.csv"))
share = data_value(data, "Share of small firms")

ge_try(v, c_e) = (for b in [(0.9, 1.3), (0.8, 1.6), (0.7, 2.2)]
                      try return ge_eq(v, c_e; bracket = b) catch e; println("ge_eq failed on $b") end
                  end; error("no GE"))
using Optim
CF_MOMENTS = [("Debt to Collateral", "avg_b2a", "b2a", true), ("Interest rate", "avg intrate", "intrate", true),
              ("Liquidation prob.", "avg_liqprob", "liqprob", true), ("CF-reliance", "CF reliance", "CF reliance", false),
              ("CF-based to total debt", "CF share", "CF share", false)]
pct(a, b) = 100 * (b / a - 1)
agg(E, r) = E.ss[r, "values"]
function cf_reform(v)
    CB = solve_eq(merge(v, (phi_af = 0.0,)))
    CR = ge_try(merge(v, (phi_af = 0.0, zeta_Rl = v.zeta_Rl / 2)), CB.c_e)
    (YL = pct(agg(CB, "YtoL"), agg(CR, "YtoL")), mass = pct(agg(CB, "totalmass"), agg(CR, "totalmass")),
     rate_pp = agg(CR, "avg intrate") - agg(CB, "avg intrate"), wage = CR.wage, gap = CR.c_e - CB.c_e)
end
lazy = load_calibration(joinpath(run_dir, "CALIB_CF.csv"))
lazy = (pc = lazy.pc, pdef_exo_l = lazy.pdef_exo_l, zeta_Rl = lazy.zeta_Rl)
t0 = time()
check = cf_reform(lazy)
println("lazy CF-only values: Y/L ", round(check.YL, digits = 3), " (DIAG cf_YL ", diag["cf_YL"], "), ", round((time() - t0) / 60, digits = 1), " min"); flush(stdout)

set_params!((phi_af = 0.0,))
cfp = (pc = (lazy.pc, 0.0, Inf), pdef_exo_l = (lazy.pdef_exo_l, 0.0, 1.0), zeta_Rl = (lazy.zeta_Rl, 0.0, Inf))
opts = Optim.Options(iterations = 300, time_limit = 3 * 3600.0, f_reltol = 1e-4, show_trace = false)
t1 = time()
Res = calibrate_fully(cfp, CF_MOMENTS, data; options = opts, n_screen = 100, n_starts = 4, n_workers = 4, max_restarts = 2,
                      model_file = MODEL_FILE, file = joinpath(run_dir, "CALIB_CF_full.csv"),
                      screen_box = (pc = (0.5 * lazy.pc, 2.0 * lazy.pc), pdef_exo_l = (0.0002, 0.05), zeta_Rl = (0.4 * lazy.zeta_Rl, 2.5 * lazy.zeta_Rl)))
fv = (pc = Res.values.pc, pdef_exo_l = Res.values.pdef_exo_l, zeta_Rl = Res.values.zeta_Rl)
println("full CF-only fit in ", round((time() - t1) / 60, digits = 1), " min: ", fv); flush(stdout)
reset_params!(); set_params!((; ext...)); set_params!(vals)
full = cf_reform(fv)
fit_err(v) = with_params(() -> begin E = solve_eq(merge(v, (phi_af = 0.0,)))
        (d2c = agg(E, "avg_b2a"), rate = agg(E, "avg intrate"), liq = agg(E, "avg_liqprob")) end, (phi_af = 0.0,))
fl, ff = fit_err(lazy), fit_err(fv)
d = (d2c = data_value(data, "Debt to Collateral"), rate = data_value(data, "Interest rate"), liq = data_value(data, "Liquidation prob."))
err(f) = sum(abs(f[k] - d[k]) / abs(d[k]) for k in keys(d))
out = DataFrame(fit = ["lazy (Search.jl)", "full (cf_refit.jl)"], pc = [lazy.pc, fv.pc], pdef_exo_l = [lazy.pdef_exo_l, fv.pdef_exo_l],
                zeta_Rl = [lazy.zeta_Rl, fv.zeta_Rl], d2c = [fl.d2c, ff.d2c], rate = [fl.rate, ff.rate], liq = [fl.liq, ff.liq],
                error = [err(fl), err(ff)], cf_YL = [check.YL, full.YL], cf_mass = [check.mass, full.mass],
                cf_rate_pp = [check.rate_pp, full.rate_pp], cf_wage = [check.wage, full.wage], full_model_YL = parse(Float64, diag["r_YL"]))
CSV.write(joinpath(SEARCH_DIR, "better_off", "$(id)_$(gtag)_cfrefit.csv"), out)
println(out)

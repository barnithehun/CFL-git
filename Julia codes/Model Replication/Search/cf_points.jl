# cf_points.jl <run id>: the CF-only zeta/2 reform at each start's final values in runs/<id>/CALIB_CF_full_runs.csv (cf_refit.jl),
# to see how much the CF-only gain depends on which of the near-equal optima is used. Writes better_off/<id>_<grid>_cfpoints.csv.
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
runs = CSV.read(joinpath(run_dir, "CALIB_CF_full_runs.csv"), DataFrame)
last_rows = combine(groupby(runs, :start), df -> df[end:end, :])
res = DataFrame(start = Int[], error = Float64[], pc = Float64[], pdef_exo_l = Float64[], zeta_Rl = Float64[], cf_YL = Float64[], cf_mass = Float64[], cf_rate_pp = Float64[], cf_wage = Float64[])
for r in eachrow(last_rows)
    v = (pc = r.pc, pdef_exo_l = r.pdef_exo_l, zeta_Rl = r.zeta_Rl)
    o = cf_reform(v)
    push!(res, (r.start, r.error, r.pc, r.pdef_exo_l, r.zeta_Rl, o.YL, o.mass, o.rate_pp, o.wage)); println(res[end, :]); flush(stdout)
end
CSV.write(joinpath(SEARCH_DIR, "better_off", "$(id)_$(gtag)_cfpoints.csv"), res)
println(res)

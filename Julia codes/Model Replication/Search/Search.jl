###########################################################################
############ CALIBRATION SEARCH: candidate settings, one run each #########
###########################################################################
# Runs the model of Baseline.jl under alternative external parameters and data targets, calibrates the internal
# parameters, computes the main results cheaply, and writes one diagnostics row per candidate, so that candidates can
# be compared in one table (tracker.csv). Nothing in Baseline.jl changes: everything is set through the overrides
# in parameters() (set_params!) and through the targets file passed to the calibration.
#
# Usage (from any folder):
#   julia +1.12.7 --project=<repo root> Search.jl run C01 C02 ...   runs the listed candidates one after the other
#   julia +1.12.7 --project=<repo root> Search.jl collect            builds tracker.csv from runs/*/DIAG.csv
# launch.ps1 <wave> [<workers>] starts several Julia processes, each with its share of the wave's candidates.
#
# candidates.csv (one row per candidate; empty cells keep the default):
#   id, wave, label      name, wave and a short description
#   calib                lazy (calibrate_lazy), lazyr (calibrate_fully from the start only: Nelder-Mead plus one restart),
#                        full (calibrate_fully: screening, 4 starts, restarts, 4 worker processes), none (keep the start values)
#   start                params (the values written in parameters()), or the id of an earlier run whose CALIB.csv gives the starts
#   stage                cheap (baseline, zeta/2 reform, decomposition, SBRA) or full (also zeta/10, perfect credit, AB only
#                        and CF only recalibrated, figures): the whole RESULTS section
#   targets              tag of the targets file: T0 = Data Replication/calib_targets.csv, else calib_targets_<tag>.csv
#                        (written by calib_targets.R with out_tag=<tag>). Its settings file sets size_by_debt (0 for
#                        size_var = assets, 1 for ltq) and d2c_cut, so model and data split firms the same way.
#   beta, delta, phi_a, kappa, rho_e, sigma_e, nul_e, size_by_debt, d2c_cut   overrides of parameters(); phi_a also sets phi_af
#   notes                free text
# Outputs, in runs/<id>/: CALIB.csv (as calibrate_lazy writes it), RESULTS.csv (tables and equilibria, as save_results),
# DIAG.csv (one row: settings, fit, small/large moments, results, the boxes) and log.txt.

using Dates, CSV, DataFrames

const SEARCH_DIR = @__DIR__
const MODEL_DIR  = normpath(joinpath(SEARCH_DIR, ".."))
# the model file: Baseline.jl, or another one in the model folder named by the environment variable SEARCH_MODEL
const MODEL_FILE = joinpath(MODEL_DIR, get(ENV, "SEARCH_MODEL", "Baseline.jl"))
const DATA_DIR   = normpath(joinpath(MODEL_DIR, "..", "..", "Data Replication"))
const RUNS_DIR   = joinpath(SEARCH_DIR, "runs")

const EXT_NAMES = (:beta, :delta, :phi_a, :kappa, :rho_e, :sigma_e, :nul_e, :size_by_debt, :d2c_cut)

read_candidates() = CSV.read(joinpath(SEARCH_DIR, "candidates.csv"), DataFrame; types = String, stripwhitespace = true)

cell(row, name) = (name in propertynames(row) && !ismissing(row[name]) && row[name] != "") ? row[name] : nothing
num(row, name) = (v = cell(row, name); v === nothing ? nothing : parse(Float64, v))

######################## collect: tracker.csv ########################
function collect_tracker()
    files = filter(isfile, [joinpath(RUNS_DIR, d, "DIAG.csv") for d in (isdir(RUNS_DIR) ? sort(readdir(RUNS_DIR)) : String[])])
    isempty(files) && (println("no DIAG.csv under $RUNS_DIR"); return)
    rows = [CSV.read(f, DataFrame; types = String) for f in files]
    cols = String[]
    for r in rows, c in names(r)
        c in cols || push!(cols, c)
    end
    t = DataFrame([c => [c in names(r) ? r[1, c] : missing for r in rows] for c in cols])
    out = joinpath(SEARCH_DIR, "tracker.csv")
    CSV.write(out, t)
    println("tracker: $(nrow(t)) candidates -> $out")
    # a short view
    short = [c for c in ["id", "wave", "label", "status", "minutes", "error", "err_size6", "s_liq", "l_liq", "s_rate", "l_rate", "s_d2c", "l_d2c",
                         "s_CFrel", "r_YL", "r_mass", "r_s_rate", "dc_M", "dc_KY", "sb_t_YL", "sb_u_YL", "sb_c_YL", "tau_interior", "score"] if c in cols]
    show(stdout, t[:, short]; allrows = true, allcols = true, eltypes = false, summary = false); println()
end

if length(ARGS) >= 1 && ARGS[1] == "collect"
    collect_tracker()
    exit()
end
length(ARGS) >= 2 && ARGS[1] == "run" || error("usage: Search.jl run <id>... | Search.jl collect")

######################## load the model ########################
ENV["GKSwstype"] = "100"                 # Plots without a display
using LinearAlgebra
BLAS.set_num_threads(1)                  # several candidates run side by side, one core each
let code = read(MODEL_FILE, String)
    marker = "########## " * "CALIBRATION"
    i = findfirst(marker, code)
    i === nothing && error("no $marker line in $MODEL_FILE")
    model = code[1:prevind(code, first(i))]
    task_local_storage(() -> include_string(Main, model, MODEL_FILE), :SOURCE_PATH, MODEL_FILE)
end
include(joinpath(MODEL_DIR, "Functions", "Results.jl"))
include(joinpath(MODEL_DIR, "Functions", "CapAlloc2.jl"))
include(joinpath(MODEL_DIR, "Functions", "Xallocation.jl"))
cd(MODEL_DIR)

# as in the calibration block of Baseline.jl
const CALIB_MOMENTS = [
    ("Debt to Collateral",     "avg_b2a",     "b2a",         true),
    ("Interest rate",          "avg intrate", "intrate",     true),
    ("Liquidation prob.",      "avg_liqprob", "liqprob",     true),
    ("CF-reliance",            "CF reliance", "CF reliance", true),
    ("CF-based to total debt", "CF share",    "CF share",    true),
]
const VALID = p -> (1 - p.phi_c_hh) + p.phi_c <= 1

# Extra items (candidates.csv column extra, ";"-separated):
#   "gap:liqprob:1", "small:CF reliance:1", "large:liqprob:1"  a size target: small minus large firms, small firms or large
#        firms' value of that sumSSsme row, added to the calibration moments (weight 1 only). Each trial then sets Fcut to
#        the data's share of small firms (Calibrate.jl, sme_at_share). They work in lazy, lazyr and full calibrations.
#   "free:sigma_v:1.0:0:5"   frees a parameter for the calibration (start, lower, upper)
#   "set:tau_tie:2"          an external override of parameters() (any number field of the model file's parameters())
#   "drop:liqprob"           the aggregate moment of that row is shown, not targeted (put a size target in its place)
#   "restarts:1"             restarts per start in a full calibration (default 3)
#   "loss:sq"                the calibration minimizes the sum of squared percentage errors (default "abs": absolute ones)
# The error reported as "error" in DIAG.csv is the five aggregate moments only (dropped ones included); "error_all" is what
# the calibration minimized.
# screening box of calibrate_fully for freed parameters (others get start x (0.5, 1.5))
const FREE_BOX = Dict(:sigma_v => (0.3, 2.0), :liq_forced => (0.02, 0.5), :sv_slope => (-0.5, 0.5))
const OPTIONS = Optim.Options(iterations = 200, time_limit = 3 * 3600.0, f_reltol = 0.002, show_trace = false)

######################## helpers ########################
# targets file of a tag, with its settings; waits for the file if calib_targets.R is still writing it
function targets_file(tag)
    file = joinpath(DATA_DIR, tag in ("", "T0") ? "calib_targets.csv" : "calib_targets_$tag.csv")
    for _ in 1:60
        isfile(file) && return file
        println("waiting for ", file); flush(stdout); sleep(30)
    end
    error("no targets file $file")
end

function targets_settings(file)
    s = replace(file, ".csv" => "_settings.csv")
    isfile(s) || return Dict{String, String}()
    t = CSV.read(s, DataFrame; types = String)
    Dict(String(r.setting) => String(r.value) for r in eachrow(t))
end

# GE wage with a wider bracket if the first one does not contain it
function ge_safe(vals, c_e; brackets = [(0.9, 1.3), (0.8, 1.6), (0.7, 2.2)])
    err = nothing
    for b in brackets
        try
            return ge_eq(vals, c_e; bracket = b)
        catch e
            err = e
            println("ge_eq failed on bracket $b: ", sprint(showerror, e))
        end
    end
    rethrow(err)
end

# share of borrowers (by mass) at an interior tau and at tau = 1
function tau_shares(E)
    S = E.SumPol; n = size(S, 1) - 1
    w = E.mu[1:n] .* (S[1:n, 4] .> 0)
    tot = sum(w)
    tau = S[1:n, 17]
    (interior = sum(w[(tau .> 0) .& (tau .< 1)]) / tot, one = sum(w[tau .== 1]) / tot)
end

# grid-edge numbers (check_grid prints its warnings itself)
function grid_flags(E)
    c = check_grid(E)
    (e_low = c.e_low, e_top_ratio = c.e_top / max(c.e_below_top, 1e-12), k_top = c.k_top, b_top = c.b_top, x_ends = max(c.x_low, c.x_high))
end

sme(E, row, col) = E.sme[row, col]
agg(E, row) = E.ss[row, "values"]
pct(b, x) = 100 * (x / b - 1)
pp(b, x) = 100 * (x - b)

######################## one candidate ########################
function run_candidate(row)
    id = row.id
    dir = joinpath(RUNS_DIR, id)
    mkpath(dir)
    t0 = now()
    D = Pair{String, Any}[]                   # the DIAG row, in order
    put!(k, v) = push!(D, k => v)
    put!("id", id); put!("wave", something(cell(row, :wave), "")); put!("label", something(cell(row, :label), ""))
    put!("status", "running"); put!("started", Dates.format(t0, "yyyy-mm-dd HH:MM"))
    calib = something(cell(row, :calib), "lazy"); stage = something(cell(row, :stage), "cheap"); start = something(cell(row, :start), "params")
    tag = something(cell(row, :targets), "T0")
    put!("calib", calib); put!("stage", stage); put!("start", start); put!("targets", tag)

    write_diag() = CSV.write(joinpath(dir, "DIAG.csv"), DataFrame([k => [v] for (k, v) in D]))

    log = open(joinpath(dir, "log.txt"), "w")
    try
        redirect_stdout(log) do
            redirect_stderr(log) do
                println("Candidate $id ($(cell(row, :label))) started ", t0); flush(log)

                #### targets and the split of small and large firms ####
                tfile = targets_file(tag)
                settings = targets_settings(tfile)
                for k in ("size_var", "cutoff_usd", "micro_usd", "d2c_cut", "end_year")
                    put!("data_" * k, get(settings, k, ""))
                end
                data_targets = read_targets(tfile)
                extra = cell(row, :extra)
                free = Pair{Symbol, Tuple{Float64, Float64, Float64}}[]     # parameters freed for the calibration
                sets = Pair{Symbol, Float64}[]                               # external overrides from extra
                moments = copy(CALIB_MOMENTS)                                # the calibration moments of this candidate
                max_restarts_full = 3                                        # restarts per start in a full calibration
                CALIB_LOSS[] = :abs                                          # error term of each moment (extra "loss:sq")
                for item in (extra === nothing ? String[] : split(extra, ";"))
                    f = split(item, ":")
                    if f[1] == "loss"
                        # "loss:sq": the calibration minimizes squared percentage errors instead of absolute ones
                        CALIB_LOSS[] = Symbol(f[2])
                        CALIB_LOSS[] in (:abs, :sq) || error("unknown loss $item")
                    elseif f[1] == "free"
                        # "free:sigma_v:0.5:0:5" (name, start, lower, upper)
                        push!(free, Symbol(f[2]) => Tuple(parse.(Float64, f[3:5])))
                    elseif f[1] == "set"
                        push!(sets, Symbol(f[2]) => parse(Float64, f[3]))
                    elseif f[1] == "restarts"
                        # "restarts:1": restarts per start in a full calibration (default 3)
                        max_restarts_full = parse(Int, f[2])
                    elseif f[1] == "agg" && f[2] == "spread"
                        # "agg:spread:1": the aggregate credit spread (sumSS avg_spread) as a target; needs a targets file
                        # with a "Credit spread" row (calib_targets.R, e.g. tag L8sp)
                        push!(moments, ("Credit spread", "avg_spread", "spread", true))
                    elseif f[1] == "drop"
                        # "drop:liqprob": the aggregate moment of that sumSS row is shown but not targeted
                        dname = Dict("liqprob" => "Liquidation prob.", "intrate" => "Interest rate", "b2a" => "Debt to Collateral",
                                     "CF reliance" => "CF-reliance", "CF share" => "CF-based to total debt", "spread" => "Credit spread")[f[2]]
                        moments = [(m, r1, r2, (m == dname && !(r1 in SIZE_ROWS)) ? false : t) for (m, r1, r2, t) in moments]
                    else
                        # "small:liqprob:1", "large:liqprob:1" or "gap:liqprob:1" (group, sumSSsme row, weight 1)
                        grp, m, w = f
                        grp in ("gap", "small", "large") || error("unknown extra item $item")
                        parse(Float64, w) == 1 || error("size targets take weight 1 only: $item")
                        dname = Dict("liqprob" => "Liquidation prob.", "intrate" => "Interest rate", "b2a" => "Debt to Collateral",
                                     "CF reliance" => "CF-reliance", "CF share" => "CF-based to total debt", "spread" => "Credit spread")[m]
                        push!(moments, (dname, String(grp), String(m), true))
                    end
                end
                put!("extra", something(extra, ""))

                #### external overrides ####
                reset_params!()
                ext = Dict{Symbol, Float64}()
                haskey(settings, "size_var") && (ext[:size_by_debt] = settings["size_var"] == "assets" ? 0.0 : 1.0)
                haskey(settings, "d2c_cut") && (ext[:d2c_cut] = parse(Float64, settings["d2c_cut"]))
                for name in EXT_NAMES
                    v = num(row, name)
                    v === nothing || (ext[name] = v)
                end
                haskey(ext, :phi_a) && !haskey(ext, :phi_af) && (ext[:phi_af] = ext[:phi_a])
                for (k, v) in sets
                    ext[k] = v; put!("x_" * String(k), v)
                end
                set_params!((; ext...))
                p = parameters()
                for name in (:beta, :delta, :phi_a, :phi_af, :kappa, :rho_e, :sigma_e, :nul_e, :DRS, :size_by_debt, :d2c_cut, :tauchen_lo)
                    put!(String(name), p[name])
                end
                put!("model_file", basename(MODEL_FILE))
                println("External overrides: ", ext); flush(log)

                #### starting values of the internal parameters ####
                # start: params, the id of an earlier run, a CALIB file, or the values themselves as v:pc=37.1;zeta_Rl=205;...
                s = start == "params" ? parameters() :
                    startswith(start, "v:") ? merge(parameters(), (; (Symbol(split(kv, "=")[1]) => parse(Float64, split(kv, "=")[2]) for kv in split(start[3:end], ";"))...)) :
                    load_calibration(isfile(start) ? start : joinpath(RUNS_DIR, start, "CALIB.csv"))
                calib_params = (pc = (s.pc, 0.0, Inf), pdef_exo_l = (s.pdef_exo_l, 0.0, 1.0), zeta_Rl = (s.zeta_Rl, 0.0, Inf),
                                phi_c_hh = (s.phi_c_hh, 0.0, 1.0), phi_c = (s.phi_c, 0.0, 1.0))
                start_vals = (pc = s.pc, pdef_exo_l = s.pdef_exo_l, zeta_Rl = s.zeta_Rl, phi_c_hh = s.phi_c_hh, phi_c = s.phi_c)
                # freed parameters (extra "free:..."): start from the start run's value if it has one, else from extra
                for (name, (st, lo, hi)) in free
                    st = (start != "params" && !startswith(start, "v:") && haskey(s, name)) ? s[name] : st
                    calib_params = merge(calib_params, NamedTuple{(name,)}(((st, lo, hi),)))
                    start_vals = merge(start_vals, NamedTuple{(name,)}((st,)))
                end
                println("Starts: ", start_vals); flush(log)

                #### calibration ####
                cfile = joinpath(dir, "CALIB.csv")
                t_cal = @elapsed CalRes =
                    calib == "lazy"  ? calibrate_lazy(calib_params, moments, data_targets; options = OPTIONS, valid = VALID, file = cfile) :
                    calib == "lazyr" ? calibrate_fully(calib_params, moments, data_targets; options = OPTIONS, valid = VALID,
                                                       starts = [start_vals], max_restarts = 1, restart_tol = 2e-3, n_workers = 0, file = cfile) :
                    calib == "full"  ? calibrate_fully(calib_params, moments, data_targets; options = OPTIONS, valid = VALID,
                                                       n_screen = 100, n_starts = 4, n_workers = 4, file = cfile, model_file = MODEL_FILE, max_restarts = max_restarts_full,
                                                       screen_box = merge((pc = (0.5 * s.pc, 2.0 * s.pc), pdef_exo_l = (0.002, 0.05), zeta_Rl = (0.4 * s.zeta_Rl, 2.5 * s.zeta_Rl),
                                                                     phi_c_hh = (0.5, 0.995), phi_c = (0.1, 0.6)),
                                                                     (; (n => FREE_BOX[n] for (n, _) in free if haskey(FREE_BOX, n))...))) :
                    calib == "none"  ? (fit = model_fit(start_vals, moments, data_targets); save_calibration(cfile, calib_params, fit);
                                        (values = fit.values, fit = fit)) :
                    error("unknown calib $calib")
                fit = CalRes.fit
                put!("calib_minutes", round(t_cal / 60, digits = 1))
                for (k, v) in pairs(fit.values)
                    put!(String(k), v)
                end
                short = Dict("Debt to Collateral" => "d2c", "Interest rate" => "rate", "Liquidation prob." => "liq",
                             "CF-reliance" => "CFrel", "CF-based to total debt" => "CFsh")
                agg_rows = [r for r in eachrow(fit.targeted) if haskey(short, r.moment)]   # the five aggregate moments
                put!("c_e", fit.c_e); put!("error", sum(abs(r.model - r.data) / abs(r.data) for r in agg_rows)); put!("error_all", fit.error)
                for r in agg_rows
                    put!("m_" * short[r.moment], r.model); put!("d_" * short[r.moment], r.data)
                end
                for r in eachrow(fit.targeted)   # credit spread target (extra "agg:spread" / "gap:spread"), if any
                    r.moment == "Credit spread" && (put!("m_spread", r.model); put!("d_spread", r.data))
                end
                for r in Iterators.filter(r -> haskey(short, r.moment), eachrow(fit.size))
                    k = short[r.moment]
                    put!("s_" * k, r.model_small); put!("l_" * k, r.model_large); put!("ds_" * k, r.data_small); put!("dl_" * k, r.data_large)
                end
                put!("share_small", fit.share.model); put!("share_data", fit.share.data)
                size6 = sum(abs(r.model_small - r.data_small) / r.data_small + abs(r.model_large - r.data_large) / r.data_large
                            for r in eachrow(fit.size) if r.moment in ("Liquidation prob.", "Interest rate", "Debt to Collateral"))
                size10 = sum(abs(r.model_small - r.data_small) / r.data_small + abs(r.model_large - r.data_large) / r.data_large for r in eachrow(fit.size) if haskey(short, r.moment))
                put!("err_size6", size6); put!("err_size10", size10)
                gap(k) = (fit.size[fit.size.moment .== k, :model_large][1] - fit.size[fit.size.moment .== k, :model_small][1]) /
                         (fit.size[fit.size.moment .== k, :data_large][1] - fit.size[fit.size.moment .== k, :data_small][1])
                for (m, k) in short
                    put!("gapratio_" * k, gap(m))
                end
                write_diag(); flush(log)

                #### results at the calibrated values ####
                reset_params!(); set_params!((; ext...)); set_params!(CalRes.values)
                BaseEq = solve_eq()
                g = grid_flags(BaseEq)
                for (k, v) in pairs(g)
                    put!("grid_" * String(k), v)
                end
                put!("base_mass", agg(BaseEq, "totalmass")); put!("base_exitshare", agg(BaseEq, "exitshare")); put!("base_meanL", agg(BaseEq, "meanL"))
                ts = tau_shares(BaseEq)
                put!("tau_interior", ts.interior); put!("tau_one", ts.one)

                ReformEq = ge_safe((zeta_Rl = parameters().zeta_Rl / 2,), BaseEq.c_e)
                put!("r_wage", ReformEq.wage); put!("r_gap_pct", 100 * (ReformEq.c_e - BaseEq.c_e) / BaseEq.c_e)
                put!("r_YL", pct(agg(BaseEq, "YtoL"), agg(ReformEq, "YtoL")))
                put!("r_mass", pct(agg(BaseEq, "totalmass"), agg(ReformEq, "totalmass")))
                put!("r_meanL", pct(agg(BaseEq, "meanL"), agg(ReformEq, "meanL")))
                put!("r_KL", pct(agg(BaseEq, "KtoL"), agg(ReformEq, "KtoL")))
                put!("r_liq_pp", pp(agg(BaseEq, "avg_liqprob"), agg(ReformEq, "avg_liqprob")))
                put!("r_d2c_pp", pp(agg(BaseEq, "avg_b2a"), agg(ReformEq, "avg_b2a")))
                put!("r_rate_pp", agg(ReformEq, "avg intrate") - agg(BaseEq, "avg intrate"))
                put!("r_CFrel_pp", pp(agg(BaseEq, "CF reliance"), agg(ReformEq, "CF reliance")))
                put!("r_CFsh_pp", pp(agg(BaseEq, "CF share"), agg(ReformEq, "CF share")))
                for (k, r) in ("liq" => "liqprob", "rate" => "intrate", "d2c" => "b2a", "CFrel" => "CF reliance", "CFsh" => "CF share")
                    put!("r_s_" * k, sme(ReformEq, r, "SME Values")); put!("r_l_" * k, sme(ReformEq, r, "LE Values"))
                end
                put!("r_share_small", sme(ReformEq, "share", "SME Values"))
                rts = tau_shares(ReformEq)
                put!("r_tau_interior", rts.interior); put!("r_tau_one", rts.one)

                Decomp = CapAlloc2(BaseEq, ReformEq)
                put!("dc_total_lp", Decomp.table.log_points[1])
                for (k, i) in ("KY" => 2, "M" => 3, "C" => 4, "Omega" => 5)
                    put!("dc_" * k, Decomp.table.share[i])
                end
                write_diag(); flush(log)

                #### SBRA ####
                SBRA = find_zeta_Rs(BaseEq; drop = 0.13)
                theta_control = parameters().phi_c * SBRA.zeta_Rs / parameters().zeta_Rl
                put!("sb_share", SBRA.zeta_Rs / parameters().zeta_Rl); put!("sb_liq_target", SBRA.target); put!("sb_liq_found", SBRA.liqprob)
                put!("theta_control", theta_control)
                TargetedEq  = ge_safe((zeta_Rs = SBRA.zeta_Rs,), BaseEq.c_e)
                UniversalEq = ge_safe((zeta_Rl = SBRA.zeta_Rs,), BaseEq.c_e)
                ControlEq   = ge_safe((zeta_Rs = SBRA.zeta_Rs, phi_c = theta_control), BaseEq.c_e)
                for (k, E) in ("t" => TargetedEq, "u" => UniversalEq, "c" => ControlEq)
                    put!("sb_$(k)_YL", pct(agg(BaseEq, "YtoL"), agg(E, "YtoL")))
                    put!("sb_$(k)_mass", pct(agg(BaseEq, "totalmass"), agg(E, "totalmass")))
                    put!("sb_$(k)_liq_pp", pp(agg(BaseEq, "avg_liqprob"), agg(E, "avg_liqprob")))
                    put!("sb_$(k)_CFsh_pp", pp(agg(BaseEq, "CF share"), agg(E, "CF share")))
                    put!("sb_$(k)_s_liq", sme(E, "liqprob", "SME Values"))
                    put!("sb_$(k)_gap_pct", 100 * (E.c_e - BaseEq.c_e) / BaseEq.c_e)
                end
                write_diag(); flush(log)

                tables = ["tab:SSvalues_firms" => size_table(["Baseline" => BaseEq, "Reform" => ReformEq]),
                          "tab:SSvalues_GE" => ge_table(BaseEq, ReformEq), "tab:decomposition" => Decomp.table,
                          "tab:percentage_changes2" => change_table(["Targeted" => (BaseEq, TargetedEq), "Universal" => (BaseEq, UniversalEq), "Control Reduc." => (BaseEq, ControlEq)]),
                          "SBRA zeta_Rs search" => SBRA.grid]
                eqs = ["Baseline" => BaseEq, "Reform" => ReformEq, "Targeted" => TargetedEq, "Universal" => UniversalEq, "Control reduction" => ControlEq]
                figs = Pair{String, Any}[]

                #### full stage: the rest of the RESULTS section ####
                if stage == "full"
                    Reform10Eq = ge_safe((zeta_Rl = parameters().zeta_Rl / 10,), BaseEq.c_e)
                    PerfEq     = ge_safe((perf_cred = 1.0,), BaseEq.c_e; brackets = [(1.0, 1.6), (1.0, 2.5)])
                    put!("r10_YL", pct(agg(BaseEq, "YtoL"), agg(Reform10Eq, "YtoL")))
                    put!("perf_YL", pct(agg(BaseEq, "YtoL"), agg(PerfEq, "YtoL")))
                    Xalloc = Xallocation(["baseline (ζ)" => BaseEq, "reform (ζ/2)" => ReformEq, "reform (ζ/10)" => Reform10Eq], "perfect credit" => PerfEq)
                    pGamPol = plot(GamPol(BaseEq.SumPol, 15, 21), GamPol(BaseEq.SumPol, 26, 21), GamPol(BaseEq.SumPol, 35, 21), layout = (1, 3), size = (900, 300), margin = 5mm)
                    write_diag(); flush(log)

                    ABCF_moments = [(moment, r1, r2, moment in ("Debt to Collateral", "Interest rate", "Liquidation prob.")) for (moment, r1, r2, _) in CALIB_MOMENTS]
                    ABCF_params = (pc = (parameters().pc, 0.0, Inf), pdef_exo_l = (parameters().pdef_exo_l, 0.0, 1.0), zeta_Rl = (parameters().zeta_Rl, 0.0, Inf))
                    ABVals = with_params(() -> calibrate_lazy(ABCF_params, ABCF_moments, data_targets; options = OPTIONS, file = joinpath(dir, "CALIB_AB.csv")).values, (phi_c = 0.0,))
                    ABBaseEq   = solve_eq(merge(ABVals, (phi_c = 0.0,)))
                    ABReformEq = ge_safe(merge(ABVals, (phi_c = 0.0, zeta_Rl = ABVals.zeta_Rl / 2)), ABBaseEq.c_e)
                    CFVals = with_params(() -> calibrate_lazy(ABCF_params, ABCF_moments, data_targets; options = OPTIONS, file = joinpath(dir, "CALIB_CF.csv")).values, (phi_af = 0.0,))
                    CFBaseEq   = solve_eq(merge(CFVals, (phi_af = 0.0,)))
                    CFReformEq = ge_safe(merge(CFVals, (phi_af = 0.0, zeta_Rl = CFVals.zeta_Rl / 2)), CFBaseEq.c_e)
                    put!("ab_YL", pct(agg(ABBaseEq, "YtoL"), agg(ABReformEq, "YtoL"))); put!("ab_mass", pct(agg(ABBaseEq, "totalmass"), agg(ABReformEq, "totalmass")))
                    put!("cf_YL", pct(agg(CFBaseEq, "YtoL"), agg(CFReformEq, "YtoL"))); put!("cf_mass", pct(agg(CFBaseEq, "totalmass"), agg(CFReformEq, "totalmass")))
                    put!("cf_rate_pp", agg(CFReformEq, "avg intrate") - agg(CFBaseEq, "avg intrate"))
                    PctTab = change_table(["AB and CF debt" => (BaseEq, ReformEq), "Only AB debt" => (ABBaseEq, ABReformEq), "Only CF debt" => (CFBaseEq, CFReformEq)])
                    push!(tables, "tab:percentage_changes" => PctTab)
                    append!(eqs, ["Reform zeta/10" => Reform10Eq, "Perfect credit" => PerfEq, "AB baseline" => ABBaseEq, "AB reform" => ABReformEq,
                                  "CF baseline" => CFBaseEq, "CF reform" => CFReformEq])
                    append!(figs, ["firm_distribution" => Xalloc.firmalloc, "capital_allocation" => Xalloc.capital, "GamPol" => pGamPol])
                end
                save_results(dir, tables, eqs, figs)

                #### the boxes (the user's list, 2026-11-06) ####
                v(k) = D[findlast(x -> x.first == k, D)].second
                boxes = Pair{String, Bool}[
                    "b1_fit"       => v("error") <= 0.30,
                    "b2_liq_s"     => abs(v("s_liq") - v("ds_liq")) <= 0.12,
                    "b2_liq_l"     => abs(v("l_liq") - v("dl_liq")) <= 0.10,
                    "b2_rategap"   => abs((v("s_rate") - v("l_rate")) - (v("ds_rate") - v("dl_rate"))) <= 1.0,
                    "b2_d2cgap"    => sign(v("l_d2c") - v("s_d2c")) == sign(v("dl_d2c") - v("ds_d2c")) && abs((v("l_d2c") - v("s_d2c")) - (v("dl_d2c") - v("ds_d2c"))) <= 0.15,
                    "b3_liq_falls" => v("r_liq_pp") < 0,
                    "b4_s_rate"    => v("r_s_rate") < v("s_rate"),
                    "b4_both"      => v("r_s_rate") < v("s_rate") && v("r_l_rate") < v("l_rate"),
                    "b5_YL"        => v("r_YL") >= 1.0,
                    "b6_mass"      => v("dc_M") >= 0.5 && v("dc_M") >= max(v("dc_KY"), v("dc_Omega")),
                    "b6_KY"        => v("dc_KY") > 0,
                    "b7_comp"      => v("dc_C") < 0,
                    "b10_target"   => v("sb_t_YL") > 0 && v("sb_t_YL") >= 0.7 * v("sb_u_YL"),
                    "b11_control"  => v("sb_c_YL") < v("sb_t_YL"),
                    "b11_negative" => v("sb_c_YL") < 0,
                    "b13_tau"      => v("tau_interior") >= 0.10,
                    "b2_CFrel_s"   => abs(v("s_CFrel") - v("ds_CFrel")) <= 0.10,
                    "b2_CFrel_l"   => abs(v("l_CFrel") - v("dl_CFrel")) <= 0.10,
                    "b5_YL12"      => v("r_YL") >= 1.2,
                ]
                for (k, b) in boxes
                    put!(k, b)
                end
                scored = ("b1_fit", "b2_liq_s", "b2_liq_l", "b2_rategap", "b2_d2cgap", "b3_liq_falls", "b4_s_rate", "b5_YL", "b6_mass", "b6_KY", "b7_comp", "b10_target", "b11_control", "b13_tau")
                put!("score", count(k -> v(k), scored))
                # required since 2026-11-07: small and large CF reliance within 0.10 of the data, and Y/L at least 1%
                put!("req_pass", v("b2_CFrel_s") && v("b2_CFrel_l") && v("b5_YL"))
                D[findfirst(x -> x.first == "status", D)] = "status" => "done"
            end
        end
    catch e
        D[findfirst(x -> x.first == "status", D)] = "status" => "error"
        put!("error_message", sprint(showerror, e)[1:min(end, 300)])
        println(log, "ERROR: ", sprint(showerror, e, catch_backtrace()))
    finally
        put!("minutes", round(Dates.value(now() - t0) / 60000, digits = 1))
        write_diag()
        println(log, "Candidate $id finished ", now(), " (", D[end].second, " min)")
        close(log)
    end
    return D
end

######################## run the listed candidates ########################
specs = read_candidates()
for id in ARGS[2:end]
    i = findfirst(==(id), specs.id)
    i === nothing && (println("no candidate $id in candidates.csv"); continue)
    println("=== $id ===")
    D = run_candidate(specs[i, :])
    println("$id: ", D[findfirst(x -> x.first == "status", D)].second, ", ", D[end].second, " min")
    flush(stdout)
end

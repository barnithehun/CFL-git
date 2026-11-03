####### CALIBRATION #######
# A calibration is written as data, in the calibration block of baseline_calib.jl:
#   data_targets:  DataFrame with columns moment, weighted, small, large - read from calib_targets.R's output (read_targets).
#                  The row "Share of small firms" sets the size cutoff Fcut.
#   calib_moments: vector of (moment, row in sumSS, row in sumSSsme, targeted); targeted = false: shown, not targeted.
#   calib_params:  NamedTuple, name = (start, lower, upper). Any number field of parameters() can be calibrated;
#                  a parameter that is not listed keeps its value from parameters().
# Every trial runs the same FirmOptim, EntryValue and sumSS as the results, with the trial values set as overrides in parameters().
# The calibration and the model fit use unrounded interest rates and spreads (rate_digits = nothing in sumSS and sumSSsme).

using Distributed   # calibrate_thoroughly runs its starts on separate Julia processes

####### Overrides in parameters() #######
# set_params!: parameters() returns vals in place of its defaults, until reset_params!()
function set_params!(vals)
    for name in keys(vals)
        (haskey(parameters(), name) && parameters()[name] isa Real) || error("parameters() has no number field :$name")
    end
    for (name, v) in pairs(vals)
        PARAM_OVERRIDE[name] = v
    end
    return parameters()
end

reset_params!() = (empty!(PARAM_OVERRIDE); parameters())

# with_params: runs f() with vals as overrides, then restores the overrides that were set before
function with_params(f, vals)
    old = copy(PARAM_OVERRIDE)
    try
        set_params!(vals)
        return f()
    finally
        empty!(PARAM_OVERRIDE)
        merge!(PARAM_OVERRIDE, old)
    end
end

####### Data targets #######
# read_targets: the table written by calib_targets.R (moment, small, large, weighted, unweighted), with its settings
function read_targets(file)
    settings = replace(file, ".csv" => "_settings.csv")
    println("Targets: ", file)
    isfile(settings) && println("  settings: ", join(["$(r.setting) = $(r.value)" for r in eachrow(CSV.read(settings, DataFrame))], ", "))
    return CSV.read(file, DataFrame)[:, [:moment, :weighted, :small, :large]]
end

# data value of a moment; col = :weighted, :small or :large
function data_value(data_targets, moment, col = :weighted)
    i = findfirst(==(moment), data_targets.moment)
    i === nothing && error("no row \"$moment\" in the data targets")
    return data_targets[i, col]
end

####### Model moments and the error #######
# error term of one moment: |model - data| / data if it is targeted, 0 if it is only shown
error_term(model, data, targeted) = targeted ? abs(model - data) / data : 0.0

# targeted moments from the sumSS table r: data, model and the error term of each moment in calib_moments
function targeted_table(r, calib_moments, data_targets)
    data  = [data_value(data_targets, moment) for (moment, _, _, _) in calib_moments]
    model = [r[row, "values"] for (_, row, _, _) in calib_moments]
    targeted = [targeted for (_, _, _, targeted) in calib_moments]
    return DataFrame(moment = first.(calib_moments), data = data, model = model,
                     error = error_term.(model, data, targeted), targeted = targeted)
end

# calib_fit: solves the model at wage with vals as overrides; returns the targeted moments, the total error and c_e
function calib_fit(vals, calib_moments, data_targets; wage = 1)
    return with_params(vals) do
        SumPol, e_chain, Fmat = FirmOptim(wage)
        c_e, f0 = EntryValue(SumPol, e_chain)
        # unrounded interest rate: sumSS's default 2 digits would make the error flat in steps
        SumRes = sumSS(SumPol, Fmat, f0; zeta_R = parameters().zeta_Rl, rate_digits = nothing)
        table = targeted_table(SumRes, calib_moments, data_targets)
        (table = table, error = sum(table.error), c_e = c_e)
    end
end

####### Objective for Nelder-Mead #######
# x holds the trial values in the order of calib_params. Outside the bounds, or when valid(parameters()) fails,
# the error is 1e10 without solving the model. ErrorFunc1 returned Inf here, but an Inf in the simplex makes
# Nelder-Mead's convergence measure NaN and stops it (the default first simplex sets phi_c_hh to 1.23).
function calib_error(x, calib_params, calib_moments, data_targets; wage = 1, valid = p -> true)
    vals = NamedTuple{keys(calib_params)}(Tuple(x))

    for (name, v) in pairs(vals)
        _, lower, upper = calib_params[name]
        lower < v < upper || return 1e10
    end
    with_params(() -> valid(parameters()), vals) || return 1e10

    println(vals)
    err = calib_fit(vals, calib_moments, data_targets; wage = wage).error
    println("error: ", err)

    return err
end

####### Size cutoff #######
# fit_Fcut: the cutoff on b' or k' (parameters().size_by_debt) that brings the model's share of small firms (the "share"
# row of sumSSsme) closest to share_target. The share only changes where Fcut passes a b' (k') value, so the current
# Fcut is kept if it gives the best split; otherwise Fcut is the value just above the split. As in sumSSsme, small firms
# are b' < Fcut among borrowers, or k' < Fcut among all firms.
function fit_Fcut(SumPol, Fmat, f0, share_target; zeta_R = parameters().zeta_Rl)
    mu, m, xpol = stat_dist(SumPol, Fmat, f0; zeta_R = zeta_R)
    n = size(SumPol, 1)
    totalmass = sum(mu)
    exitshare = (transpose(mu) * xpol) / totalmass
    by_debt = parameters().size_by_debt == 1
    kp = SumPol[1:n-1, by_debt ? 4 : 3]
    pop = by_debt ? SumPol[1:n-1, 4] .!= 0 : trues(n-1)
    denom = by_debt ? transpose(mu[1:n-1]) * pop / totalmass : 1 - exitshare
    share(F) = sum(mu[1:n-1][pop .& (kp .< F)]) / totalmass / denom

    ks = sort(unique(kp[pop]))
    cuts = [ks[2:end]; 2 * ks[end]]           # cut j: small = {k' <= ks[j]}
    gaps = [abs(share(F) - share_target) for F in cuts]
    Fcut = parameters().Fcut
    abs(share(Fcut) - share_target) <= minimum(gaps) && return Fcut
    return cuts[argmin(gaps)]
end

####### Model fit #######
# model_fit: solves the model at wage with vals as overrides, sets Fcut to match the share of small firms, and returns
# the total error, the share of small firms, the targeted moments and the small and large firm moments.
# vals = NamedTuple() gives the fit at the values in parameters().
function model_fit(vals, calib_moments, data_targets; wage = 1)
    return with_params(vals) do
        p = parameters()
        p.pdef_exo_s == p.pdef_exo_l || @warn "pdef_exo_s differs from pdef_exo_l, so Fcut also enters FirmOptim; fit_Fcut does not account for that"
        SumPol, e_chain, Fmat = FirmOptim(wage)
        c_e, f0 = EntryValue(SumPol, e_chain)
        targeted = targeted_table(sumSS(SumPol, Fmat, f0; zeta_R = p.zeta_Rl, rate_digits = nothing), calib_moments, data_targets)

        # the cutoff only moves if the targets have a share of small firms
        share_data = "Share of small firms" in data_targets.moment ? data_value(data_targets, "Share of small firms") : missing
        Fcut = ismissing(share_data) ? p.Fcut : fit_Fcut(SumPol, Fmat, f0, share_data; zeta_R = p.zeta_Rl)
        sme = with_params(() -> sumSSsme(SumPol, Fmat, f0; zeta_R = p.zeta_Rl, rate_digits = nothing), (Fcut = Fcut,))

        size_table = DataFrame(moment = first.(calib_moments),
                               data_small  = [data_value(data_targets, moment, :small) for (moment, _, _, _) in calib_moments],
                               model_small = [sme[row, "SME Values"] for (_, _, row, _) in calib_moments],
                               data_large  = [data_value(data_targets, moment, :large) for (moment, _, _, _) in calib_moments],
                               model_large = [sme[row, "LE Values"] for (_, _, row, _) in calib_moments])

        (values = merge(vals, (Fcut = Fcut,)), error = sum(targeted.error), c_e = c_e, wage = wage, Fcut = Fcut,
         share = (data = share_data, model = sme["share", "SME Values"]), targeted = targeted, size = size_table)
    end
end

function print_fit(fit)
    println("\nTotal error (minimized by calibrate): ", round(fit.error, digits = 5), "     c_e = ", round(fit.c_e, digits = 2), " at wage ", fit.wage)
    println("Share of small firms: data ", round(fit.share.data, digits = 3), ", model ", round(fit.share.model, digits = 3),
            "  (small: ", parameters().size_by_debt == 1 ? "b' < Fcut among borrowers" : "k' < Fcut", ", Fcut = ", round(fit.Fcut, digits = 2), ")")
    println("\nTargeted moments, all firms (targeted = false: shown, not targeted)")
    show(stdout, fit.targeted; allrows = true, allcols = true, eltypes = false); println()
    println("\nSmall and large firms (not targeted)")
    show(stdout, fit.size; allrows = true, allcols = true, eltypes = false); println()
end

####### Saving #######
# save_calibration: one CSV, by section - internally calibrated parameters, Fcut, externally set parameters, grid,
# the moments, the share of small firms and the fit
function save_calibration(file, calib_params, fit)
    internal = [keys(calib_params)...]
    p, g = with_params(() -> (parameters(), gridsize()), fit.values)

    row(section, name, value; data = missing, start = missing, error = missing) =
        (section = section, name = String(name), value = value, data = data, start = start, error = error)
    rows = NamedTuple[]
    append!(rows, [row("internal parameter", name, p[name]; start = calib_params[name][1]) for name in internal])
    push!(rows, row("Fcut, set after the fit", :Fcut, fit.Fcut; data = fit.share.data))
    for name in keys(p)
        name in internal || name == :Fcut || push!(rows, row("external parameter", name, name == :tau_vec ? string(p[name]) : p[name]))
    end
    append!(rows, [row("grid", name, g[name]) for name in keys(g)])
    append!(rows, [row(r.targeted ? "targeted moment" : "moment, not targeted", r.moment, r.model; data = r.data, error = r.error)
                   for r in eachrow(fit.targeted)])
    append!(rows, [row("small firms", r.moment, r.model_small; data = r.data_small) for r in eachrow(fit.size)])
    append!(rows, [row("large firms", r.moment, r.model_large; data = r.data_large) for r in eachrow(fit.size)])
    push!(rows, row("share of small firms", "Share of small firms", fit.share.model; data = fit.share.data))
    append!(rows, [row("fit", "total error", fit.error), row("fit", "c_e", fit.c_e), row("fit", "wage", fit.wage)])

    mkpath(dirname(file))
    CSV.write(file, DataFrame(rows))
    println("Saved to ", file)
end

####### CALIBRATE #######
# Runs Nelder-Mead from the starts in calib_params, then sets Fcut to match the share of small firms. Prints the
# lines to paste into parameters(), saves the calibration to file and returns the calibrated values (Fcut included),
# the model fit and the Optim result.
function calibrate(calib_params, calib_moments, data_targets; options = Optim.Options(), wage = 1, valid = p -> true,
                   file = joinpath("calibration_outputs", "CALIB_$(Dates.format(now(), "yyyy-mm-dd_HHMM")).csv"))

    x0 = [start for (start, _, _) in values(calib_params)]
    OptRes = optimize(x -> calib_error(x, calib_params, calib_moments, data_targets; wage = wage, valid = valid), x0, NelderMead(), options)
    println(OptRes)

    vals = NamedTuple{keys(calib_params)}(Tuple(Optim.minimizer(OptRes)))
    fit = model_fit(vals, calib_moments, data_targets; wage = wage)

    println("\nTo keep these values, replace these lines in parameters():")
    for (name, v) in pairs(fit.values)
        println("    $name::Float64 = par(:$name, $v)")
    end
    save_calibration(file, calib_params, fit)

    return (values = fit.values, fit = fit, optim = OptRes)
end

####### CALIBRATE THOROUGHLY #######
# calibrate_thoroughly: calibrate with a screening step, several starts run in parallel, restarts, and the parameters
# searched in transformed form. The error has many local minima (moments jump where a policy switches grid points),
# and in the tests on 2026-10-31 the result depended mostly on where Nelder-Mead started.
#  1. Screening: the error at n_screen points spread evenly over screen_box (a Halton sequence) and at the start in
#     calib_params, one model solve each. The best n_starts points become the starts. screen_box: name = (low, high);
#     a parameter left out gets start × (0.5, 1.5), cut to its bounds. Points outside the bounds or failing valid are
#     dropped. starts = [CalRes.values, ...] skips the screening (fields not in calib_params, such as Fcut, are ignored).
#  2. From each start, Nelder-Mead runs, then restarts from where it stopped with a new first simplex, until a restart
#     improves the error by less than restart_tol, or after max_restarts. The best result over all starts is kept.
#  3. transform = true: Nelder-Mead searches z = logit of a parameter with two finite bounds and z = log(p - lower) of
#     one with only a lower bound, so no trial leaves the bounds; the first simplex moves each z by simplex_step
#     (0.4: ×1.5 for a log, 0.80 -> 0.86 for a logit). transform = false searches the parameters themselves, with
#     Optim's default first simplex (+50% on each), as calibrate.
#  n_workers: separate Julia processes the screening and the starts run on (0 = this process only); the overrides in
#     parameters() are global, so threads would overwrite each other's. Each worker loads the model part of model_file
#     (everything above "########## CALIBRATION"), so changes to gridsize() or the model must be saved there first.
#     Overrides set with set_params! before the call are passed on to the workers.
# Saves the same CSV as calibrate (CALIB_<date>_thorough.csv; its start column shows the best run's start),
# <file>_runs.csv with every run and <file>_screen.csv with every screened point.

# Halton sequence: point i in the dimension with prime base, in (0, 1)
function halton(i, base)
    f, r = 1.0, 0.0
    while i > 0
        f /= base
        r += f * (i % base)
        i ÷= base
    end
    return r
end
const HALTON_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)

# transformed parameter z: logit((p - lower) / (upper - lower)) with two finite bounds, log(p - lower) or log(upper - p)
# with one, p itself with none
to_z(p, lower, upper) = isfinite(lower) && isfinite(upper) ? log((p - lower) / (upper - p)) :
                        isfinite(lower) ? log(p - lower) : isfinite(upper) ? log(upper - p) : p
from_z(z, lower, upper) = isfinite(lower) && isfinite(upper) ? lower + (upper - lower) / (1 + exp(-z)) :
                          isfinite(lower) ? lower + exp(z) : isfinite(upper) ? upper - exp(z) : z

# first simplex of Nelder-Mead: the start, and the start with each coordinate moved by step
struct StepSimplexer <: Optim.Simplexer
    step::Float64
end
Optim.simplexer(S::StepSimplexer, x0) = [[copy(x0)]; [x0 .+ S.step .* (eachindex(x0) .== j) for j in eachindex(x0)]]

# default screening box: start × (0.5, 1.5), cut to the bounds
default_box(calib_params) = map(((start, lower, upper),) -> (max(min(0.5start, 1.5start), lower), min(max(0.5start, 1.5start), upper)), calib_params)

# error at the parameter values vals (a NamedTuple); 1e10 outside the bounds, when valid fails or when the error is NaN.
# base: the overrides set before calibrate_thoroughly was called.
function thorough_error(vals, calib_params, calib_moments, data_targets, wage, valid, base)
    all(calib_params[n][2] < v < calib_params[n][3] for (n, v) in pairs(vals)) || return 1e10
    return with_params(base) do
        with_params(() -> valid(parameters()), vals) || return 1e10
        err = calib_fit(vals, calib_moments, data_targets; wage = wage).error
        isnan(err) ? 1e10 : err
    end
end

# thorough_start: Nelder-Mead with restarts from start x (a NamedTuple), with the settings s; runs on a worker.
# Returns every run, and the best values, error and Optim result.
function thorough_start(i, x, s)
    pnames = keys(s.calib_params)
    lower = [s.calib_params[n][2] for n in pnames]
    upper = [s.calib_params[n][3] for n in pnames]
    to_vals(z) = NamedTuple{pnames}(Tuple(s.transform ? from_z.(z, lower, upper) : z))
    f(z) = thorough_error(to_vals(z), s.calib_params, s.calib_moments, s.data_targets, s.wage, s.valid, s.base)
    method = s.transform ? NelderMead(initial_simplex = StepSimplexer(s.simplex_step)) : NelderMead()
    z0 = s.transform ? to_z.(collect(x), lower, upper) : collect(x)

    runs = NamedTuple[]
    OptRes = nothing
    for r in 0:s.max_restarts
        t = @elapsed res = optimize(f, r == 0 ? z0 : Optim.minimizer(OptRes), method, s.options)
        gain = r == 0 ? missing : Optim.minimum(OptRes) - Optim.minimum(res)
        push!(runs, merge((start = i, restart = r, error = Optim.minimum(res), gain = gain,
                           trials = Optim.f_calls(res), minutes = t / 60), to_vals(Optim.minimizer(res))))
        println("start $i, ", r == 0 ? "first run" : "restart $r", ": error ", Optim.minimum(res),
                r == 0 ? "" : ", gain $gain", ", ", Optim.f_calls(res), " trials, ", round(t / 60, digits = 1), " min")
        (r == 0 || gain > 0) && (OptRes = res)
        r > 0 && gain < s.restart_tol && break
    end
    return (runs = runs, values = to_vals(Optim.minimizer(OptRes)), error = Optim.minimum(OptRes), start = x, optim = OptRes)
end

# start_workers: n Julia processes with the model part of model_file loaded
function start_workers(n, model_file)
    code = read(model_file, String)
    i = findfirst("########## CALIBRATION", code)
    i === nothing && error("no \"########## CALIBRATION\" line in $model_file")
    model = code[1:prevind(code, first(i))]
    ws = addprocs(n; exeflags = "--project=$(Base.active_project())")
    println("Loading the model from $model_file on $n workers")
    @everywhere ws task_local_storage(() -> include_string(Main, $model, $model_file), :SOURCE_PATH, $model_file)
    return ws
end

function calibrate_thoroughly(calib_params, calib_moments, data_targets; options = Optim.Options(), wage = 1, valid = p -> true,
                              screen_box = NamedTuple(), n_screen = 100, n_starts = 4, starts = nothing,
                              restart_tol = 1e-6, max_restarts = 3, transform = true, simplex_step = 0.4,
                              n_workers = n_starts, model_file = normpath(joinpath(@__DIR__, "..", "..", "baseline_calib.jl")),
                              file = joinpath("calibration_outputs", "CALIB_$(Dates.format(now(), "yyyy-mm-dd_HHMM"))_thorough.csv"))

    pnames = keys(calib_params)
    base = (; PARAM_OVERRIDE...)
    inside(vals) = all(calib_params[n][2] < v < calib_params[n][3] for (n, v) in pairs(vals)) &&
                   with_params(() -> valid(parameters()), vals)
    s = (calib_params = calib_params, calib_moments = calib_moments, data_targets = data_targets, wage = wage, valid = valid,
         base = base, options = options, transform = transform, simplex_step = simplex_step,
         restart_tol = restart_tol, max_restarts = max_restarts)

    ws = n_workers > 0 ? start_workers(n_workers, model_file) : Int[]
    run_map(f, xs) = isempty(ws) ? map(f, xs) : pmap(f, WorkerPool(ws), xs)
    screen, results = nothing, nothing
    try
        if starts === nothing
            box = merge(default_box(calib_params), screen_box)
            points = [NamedTuple{pnames}(Tuple(box[n][1] + halton(j, HALTON_BASES[d]) * (box[n][2] - box[n][1])
                                               for (d, n) in enumerate(pnames))) for j in 1:n_screen]
            pushfirst!(points, NamedTuple{pnames}(Tuple(start for (start, _, _) in values(calib_params))))
            points = filter(inside, points)
            println("\n####### Screening $(length(points)) points #######")
            t = @elapsed errs = run_map(x -> thorough_error(x, calib_params, calib_moments, data_targets, wage, valid, base), points)
            screen = sort(DataFrame([merge((point = j, error = e), x) for (j, (x, e)) in enumerate(zip(points, errs))]), :error)
            println("Screening took ", round(t / 60, digits = 1), " min; the best points:")
            show(stdout, first(screen, n_starts); allcols = true, eltypes = false); println()
            starts = [NamedTuple{pnames}(Tuple(r[n] for n in pnames)) for r in eachrow(first(screen, n_starts))]
        else
            starts = [NamedTuple{pnames}(Tuple(Float64(x[n]) for n in pnames)) for x in starts]
            for (i, x) in enumerate(starts)
                inside(x) || @warn "start $i is outside the bounds or fails valid, skipped" start = x
            end
            starts = filter(inside, starts)
        end
        isempty(starts) && error("no start lies inside the bounds and passes valid")

        println("\n####### Nelder-Mead from $(length(starts)) starts #######")
        results = run_map(((i, x),) -> thorough_start(i, x, s), collect(enumerate(starts)))
    finally
        isempty(ws) || rmprocs(ws)
    end
    best = results[argmin([r.error for r in results])]

    runs = DataFrame(reduce(vcat, [r.runs for r in results]))
    runs.best = (1:nrow(runs)) .== argmin(runs.error)
    println("\nAll runs (if the starts end at different values with similar errors, the moments do not pin down the parameters)")
    show(stdout, runs; allrows = true, allcols = true, eltypes = false); println()

    # from here on as in calibrate
    vals = best.values
    fit = model_fit(vals, calib_moments, data_targets; wage = wage)

    println("\nTo keep these values, replace these lines in parameters():")
    for (name, v) in pairs(fit.values)
        println("    $name::Float64 = par(:$name, $v)")
    end
    # the start column of the CSV shows the start of the best run
    save_calibration(file, NamedTuple{pnames}(Tuple((best.start[n], calib_params[n][2], calib_params[n][3]) for n in pnames)), fit)
    runs_file = replace(file, ".csv" => "_runs.csv")
    CSV.write(runs_file, runs)
    println("Saved the runs to ", runs_file)
    if screen !== nothing
        screen_file = replace(file, ".csv" => "_screen.csv")
        CSV.write(screen_file, screen)
        println("Saved the screening to ", screen_file)
    end

    return (values = fit.values, fit = fit, optim = best.optim, runs = runs, screen = screen)
end

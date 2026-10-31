####### CALIBRATION #######
# A calibration is written as data, in the calibration block of baseline_calib.jl:
#   data_targets:  DataFrame with columns moment, weighted, small, large - read from calib_targets.R's output (read_targets)
#                  or written by hand. The row "Share of small firms" sets the size cutoff Fcut.
#   calib_moments: vector of (moment, row in sumSS, row in sumSSsme, error form); error form :none = shown, not targeted.
#   calib_params:  NamedTuple, name = (start, lower, upper). Any number field of parameters() can be calibrated;
#                  a parameter that is not listed keeps its value from parameters().
# Every trial runs the same FirmOptim, EntryValue and sumSS as the results, with the trial values set as overrides in parameters().
# The calibration and the model fit use unrounded interest rates and spreads (rate_digits = nothing in sumSS and sumSSsme).

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
# error term of one moment:  :rel     |model - data| / data
#                            :rel_1m  |model - data| / (1 - data), the same as comparing 1 - model with 1 - data
#                                     (ErrorFunc1 compared the AB share instead of the CF share, for greater weight)
#                            :none    0, the moment is shown but not targeted
function error_term(model, data, form)
    form == :rel    && return abs(model - data) / data
    form == :rel_1m && return abs(model - data) / (1 - data)
    form == :none   && return 0.0
    error("unknown error form :$form")
end

# targeted moments from the sumSS table r: data, model and the error term of each moment in calib_moments
function targeted_table(r, calib_moments, data_targets)
    data  = [data_value(data_targets, moment) for (moment, _, _, _) in calib_moments]
    model = [r[row, "values"] for (_, row, _, _) in calib_moments]
    form  = [form for (_, _, _, form) in calib_moments]
    return DataFrame(moment = first.(calib_moments), data = data, model = model,
                     error = error_term.(model, data, form), form = string.(form))
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
# fit_Fcut: the cutoff on k' that brings the model's share of small firms (the "share" row of sumSSsme) closest to
# share_target. The share only changes where Fcut passes a k' value, so the current Fcut is kept if it gives the best
# split; otherwise Fcut is the k' value just above the split (small firms: k' < Fcut, as in sumSSsme).
function fit_Fcut(SumPol, Fmat, f0, share_target; zeta_R = parameters().zeta_Rl)
    mu, m, xpol = stat_dist(SumPol, Fmat, f0; zeta_R = zeta_R)
    n = size(SumPol, 1)
    totalmass = sum(mu)
    exitshare = (transpose(mu) * xpol) / totalmass
    kp = SumPol[1:n-1, 3]
    share(F) = sum(mu[1:n-1][kp .< F]) / totalmass / (1 - exitshare)

    ks = sort(unique(kp))
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
            "  (small: k' < Fcut = ", round(fit.Fcut, digits = 2), ")")
    println("\nTargeted moments, all firms (form :none = not targeted)")
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
    append!(rows, [row("targeted moment ($(r.form))", r.moment, r.model; data = r.data, error = r.error) for r in eachrow(fit.targeted)])
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
# default_starts: the start in calib_params, and two variants that move each start by ±pct with alternating signs
# (+,-,+,... and -,+,-,...). A start of 0 does not move.
function default_starts(calib_params, pct)
    x0 = [start for (start, _, _) in values(calib_params)]
    variant(s) = NamedTuple{keys(calib_params)}(Tuple(x0 .* (1 .+ s .* pct)))
    signs = [isodd(i) ? 1 : -1 for i in eachindex(x0)]
    pct == 0 && return [variant(0)]
    return [variant(0), variant(signs), variant(-signs)]
end

# calibrate_thoroughly: calibrate with restarts and several starts. Nelder-Mead often stops before the minimum, when its
# simplex shrinks in a rugged part of the error. A restart runs it again from where it stopped, with a new first simplex
# (Optim's default: +50% on each parameter). Restarts stop when one improves the error by less than restart_tol, or after
# max_restarts. A start that ends in a wrong local minimum stays there, so the loop runs from each start and keeps the best.
#   starts: NamedTuples with the names of calib_params (other fields, such as Fcut, are ignored). The default is
#           default_starts(calib_params, variant_pct); starts = [CalRes.values] polishes an earlier calibrate result.
# The first run from the start in calib_params is the same as calibrate's, so the result is never worse.
# Saves the same CSV as calibrate, and <file>_runs.csv with every run: start, restart, error, gain, trials, minutes, values.
function calibrate_thoroughly(calib_params, calib_moments, data_targets; options = Optim.Options(), wage = 1, valid = p -> true,
                              variant_pct = 0.10, starts = default_starts(calib_params, variant_pct),
                              restart_tol = 1e-6, max_restarts = 5,
                              file = joinpath("calibration_outputs", "CALIB_$(Dates.format(now(), "yyyy-mm-dd_HHMM"))_thorough.csv"))

    pnames = keys(calib_params)
    f(x) = calib_error(x, calib_params, calib_moments, data_targets; wage = wage, valid = valid)
    inside(vals) = all(calib_params[n][2] < v < calib_params[n][3] for (n, v) in pairs(vals)) &&
                   with_params(() -> valid(parameters()), vals)

    runs = NamedTuple[]
    best = nothing
    for (i, s) in enumerate(starts)
        x = NamedTuple{pnames}(Tuple(Float64(s[n]) for n in pnames))
        if !inside(x)
            @warn "start $i is outside the bounds or fails valid, skipped" start = x
            continue
        end

        OptRes = nothing
        for r in 0:max_restarts
            println("\n####### start $i of $(length(starts)), ", r == 0 ? "first run" : "restart $r", " #######")
            t = @elapsed res = optimize(f, r == 0 ? collect(x) : Optim.minimizer(OptRes), NelderMead(), options)
            gain = r == 0 ? missing : Optim.minimum(OptRes) - Optim.minimum(res)
            push!(runs, merge((start = i, restart = r, error = Optim.minimum(res), gain = gain,
                               trials = Optim.f_calls(res), minutes = t / 60),
                              NamedTuple{pnames}(Tuple(Optim.minimizer(res)))))
            println("error ", Optim.minimum(res), r == 0 ? "" : ", gain $gain", ", ", Optim.f_calls(res), " trials, ",
                    round(t / 60, digits = 1), " min")
            (r == 0 || gain > 0) && (OptRes = res)
            r > 0 && gain < restart_tol && break
        end

        if best === nothing || Optim.minimum(OptRes) < Optim.minimum(best.optim)
            best = (start = x, optim = OptRes)
        end
    end
    best === nothing && error("no start lies inside the bounds and passes valid")

    runs = DataFrame(runs)
    runs.best = (1:nrow(runs)) .== argmin(runs.error)
    println("\nAll runs (if the starts end at different values with similar errors, the moments do not pin down the parameters)")
    show(stdout, runs; allrows = true, allcols = true, eltypes = false); println()

    # from here on as in calibrate
    vals = NamedTuple{pnames}(Tuple(Optim.minimizer(best.optim)))
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

    return (values = fit.values, fit = fit, optim = best.optim, runs = runs)
end

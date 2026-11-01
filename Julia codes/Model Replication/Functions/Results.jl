####### RESULTS #######
# Helpers for the RESULTS section of baseline_calib.jl and for MiscResults.jl.
# An experiment is a set of overrides in parameters() (with_params, Calibrate.jl) plus a wage, e.g. (zeta_Rl = 1376.6,)
# for the reform. Overrides add to those already set (set_params!), so the results run at the calibrated values.
# Interest rates and spreads are unrounded (rate_digits = nothing in sumSS and sumSSsme); tables round only when printed.

####### Saved calibrations #######
# load_calibration: the internally calibrated parameters and Fcut from a CSV written by calibrate or calibrate_thoroughly
function load_calibration(file)
    t = CSV.read(file, DataFrame; types = Dict(:value => String))
    rows = filter(r -> r.section == "internal parameter" || startswith(r.section, "Fcut"), t)
    println("Calibration: ", file)
    return (; (Symbol(r.name) => parse(Float64, r.value) for r in eachrow(rows))...)
end

####### Equilibrium #######
# This is a flexible setup of my model model solution
#  solve_eq: solves the model at wage with vals, where vals is a NamedTuple of overrides in parameters()
#  the last block summarizes the results in a NamedTuple 
function solve_eq(vals = NamedTuple(); wage = 1.0)
    return with_params(vals) do
        p = parameters()

        SumPol, e_chain, Fmat = FirmOptim(wage)
        c_e, f0 = EntryValue(SumPol, e_chain)
        mu, m, xpol = stat_dist(SumPol, Fmat, f0; zeta_R = p.zeta_Rl)

        (values = (; PARAM_OVERRIDE...), wage = wage, params = p, gridsize = gridsize(),
         SumPol = SumPol, e_chain = e_chain, Fmat = Fmat, f0 = f0, c_e = c_e, mu = mu, m = m, xpol = xpol,
         ss  = sumSS(SumPol, Fmat, f0; zeta_R = p.zeta_Rl, rate_digits = nothing),
         sme = sumSSsme(SumPol, Fmat, f0; zeta_R = p.zeta_Rl, rate_digits = nothing))
    end
end

# Finding the GE wage is a bisection search for the wage where the entry value equals c_e_target,
# With vals as overrides in parameters().
function ge_eq(vals, c_e_target; bracket = (0.9, 1.3), wage_tol = 1e-6, gap_tol = 1e-3)
    
    # gap(w): entry value at wage w (with vals) minus baseline c_e; bisection finds w where gap = 0
    gap(w) = with_params(() -> FindWage(w, parameters().phi_c, parameters().zeta_Rl), vals) - c_e_target

    # bisection to find the wage where gap = 0
    lo, hi = bracket
    g_lo, g_hi = gap(lo), gap(hi)
    sign(g_lo) != sign(g_hi) ||
        error("entry value - c_e has the same sign at both ends of the wage bracket $bracket ($g_lo, $g_hi); widen the bracket")
    t = @elapsed while hi - lo > wage_tol && g_lo != 0 && g_hi != 0
        mid = (lo + hi) / 2
        g_mid = gap(mid)
        if sign(g_mid) == sign(g_lo)
            lo, g_lo = mid, g_mid
        else
            hi, g_hi = mid, g_mid
        end
    end
    wage = abs(g_lo) <= abs(g_hi) ? lo : hi
    eq = solve_eq(vals; wage = wage)
    
    # prints the GE wage, the solve time and the remaining gap (entry value - c_e)
    # warns if |gap| > 0.1%
    println("GE wage ", round(wage, digits = 6), " (", round(t, digits = 1), " s): entry value - c_e = ", eq.c_e - c_e_target)
    abs(eq.c_e - c_e_target) > gap_tol * c_e_target &&
        @warn "No wage gives free entry exactly: entry value - c_e jumps from $g_lo at w = $lo to $g_hi at w = $hi. Kept w = $wage, the end closer to free entry."
    return eq
end

####### SBRA: targeted cost cut #######
# find_zeta_Rs: the fixed reorganization cost of small firms (k' < Fcut) that lowers small firms' liquidation probability
# by drop (main5: about 13 pp) from the baseline B, at B's wage. It is a grid search over shares of B's zeta_Rl, as in
# SmallCosts.jl, which found that bisection fails here (the probability is not monotone in zeta_Rs).
function find_zeta_Rs(B; drop = 0.13, shares = 0.40:0.01:1.0)
    zeta = B.params.zeta_Rl
    target = B.sme["liqprob", "SME Values"] - drop
    liq = map(shares) do s
        with_params(merge(B.values, (zeta_Rs = s * zeta,))) do
            SumPol, e_chain, Fmat = FirmOptim(B.wage)
            c_e, f0 = EntryValue(SumPol, e_chain)
            sumSSsme(SumPol, Fmat, f0; zeta_R = zeta, rate_digits = nothing)["liqprob", "SME Values"]
        end
    end
    i = argmin(abs.(liq .- target))
    println("zeta_Rs = ", round(shares[i] * zeta, digits = 2), " (", shares[i], " of zeta_Rl): small firms' liquidation probability ",
            round(liq[i], digits = 4), ", target ", round(target, digits = 4))
    return (zeta_Rs = shares[i] * zeta, liqprob = liq[i], target = target,
            grid = DataFrame(zeta_Rs = shares .* zeta, liqprob_small = liq))
end

####### Summary Tables #######
# size_table: small and large firm moments in each equilibrium (tab:SSvalues_firms); eqs = ["name" => equilibrium, ...]
function size_table(eqs)
    rows = ["Pr. Liquidation" => "liqprob", "Debt/Collateral" => "b2a", "Interest Rate" => "intrate",
            "CF Reliance" => "CF reliance", "CF-based/Total Debt" => "CF share"]
    t = DataFrame(moment = first.(rows))
    for (name, eq) in eqs
        t[!, "$name small"] = [eq.sme[row, "SME Values"] for row in last.(rows)]
        t[!, "$name large"] = [eq.sme[row, "LE Values"] for row in last.(rows)]
    end
    return t
end

# ge_table: aggregate moments in the baseline B and the reform R, and the % change (tab:SSvalues_GE)
function ge_table(B, R)
    rows = ["Total Firm Mass" => "totalmass", "Mean Employees" => "meanL",
            "Capital Intensity (K/L)" => "KtoL", "Productivity (Y/L)" => "YtoL"]
    base   = [B.ss[row, "values"] for row in last.(rows)]
    reform = [R.ss[row, "values"] for row in last.(rows)]
    return DataFrame(moment = first.(rows), Baseline = base, Reform = reform, change_pct = 100 .* (reform ./ base .- 1))
end

# change_table: % and percentage-point changes from a baseline, one column per experiment
# (tab:percentage_changes and tab:percentage_changes2); cols = ["name" => (baseline, experiment), ...]
function change_table(cols)
    rows = [("Total firm mass (%)", :pct, "totalmass"), ("Mean employees (%)", :pct, "meanL"),
            ("Productivity (Y/L) (%)", :pct, "YtoL"), ("Capital Int. (K/L) (%)", :pct, "KtoL"), ("Wage (%)", :pct, "wage"),
            ("Pr. Liquidation (pp)", :pp, "avg_liqprob"), ("Debt/Collateral (pp)", :pp, "avg_b2a"),
            ("Interest rate (pp)", :diff, "avg intrate"), ("CF reliance (pp)", :pp, "CF reliance"), ("CF share (pp)", :pp, "CF share")]
    value(eq, row) = row == "wage" ? eq.wage : eq.ss[row, "values"]
    change(b, x, kind) = kind == :pct ? 100 * (x / b - 1) : kind == :pp ? 100 * (x - b) : x - b   # the rate is already in %
    t = DataFrame(moment = [r[1] for r in rows])
    for (name, (B, X)) in cols
        t[!, name] = [change(value(B, row), value(X, row), kind) for (_, kind, row) in rows]
    end
    return t
end

function print_table(title, t; digits = 3)
    println("\n", title)
    show(stdout, transform(t, names(t, Number) .=> (c -> round.(c, digits = digits)); renamecols = false);
         allrows = true, allcols = true, eltypes = false, summary = false)
    println()
end

####### Saving #######
# save_results: dir/RESULTS.csv with one row per table cell (section = table, row, column, value) and per equilibrium
# (section "equilibrium", row = its name, column = an override, the wage or c_e), and the figures as dir/<name>.png.
# load_eqs reads the equilibria back, so MiscResults.jl can solve them again.
function save_results(dir, tables, eqs, figs)
    rows = NamedTuple[]
    for (name, t) in tables, r in eachrow(t), col in names(t)[2:end]
        push!(rows, (section = name, row = string(r[1]), column = col, value = Float64(r[col])))
    end
    for (name, eq) in eqs
        for (k, v) in pairs(eq.values)
            push!(rows, (section = "equilibrium", row = name, column = String(k), value = v))
        end
        push!(rows, (section = "equilibrium", row = name, column = "wage", value = eq.wage))
        push!(rows, (section = "equilibrium", row = name, column = "c_e", value = eq.c_e))
    end
    mkpath(dir)
    CSV.write(joinpath(dir, "RESULTS.csv"), DataFrame(rows))
    for (name, fig) in figs
        savefig(fig, joinpath(dir, "$name.png"))
    end
    println("Saved to ", dir)
end

# load_eqs: the equilibria saved by save_results, as name => (values, wage, c_e)
function load_eqs(file)
    t = filter(r -> r.section == "equilibrium", CSV.read(file, DataFrame))
    eqs = Dict{String, NamedTuple}()
    for name in unique(t.row)
        r = t[t.row .== name, :]
        over = r[.!in.(r.column, Ref(("wage", "c_e"))), :]
        eqs[name] = (values = (; (Symbol(c) => v for (c, v) in zip(over.column, over.value))...),
                     wage = only(r.value[r.column .== "wage"]), c_e = only(r.value[r.column .== "c_e"]))
    end
    return eqs
end

# newest_results: dir/<date>/RESULTS.csv of the latest run (folder names are dates, so the last in order is the latest)
function newest_results(dir = "results_outputs")
    runs = sort(filter(d -> isfile(joinpath(dir, d, "RESULTS.csv")), readdir(dir)))
    isempty(runs) && error("no RESULTS.csv in $dir; run the RESULTS section of baseline_calib.jl first")
    return joinpath(dir, runs[end], "RESULTS.csv")
end

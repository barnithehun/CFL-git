############ Steady-state consumption (welfare): baseline, zeta/2 and zeta/3 reforms #########
# Loads a model file up to its CALIBRATION block (default Baseline_noise60.jl, ENV AGG_MODEL to change), solves the
# baseline and the two reforms at their GE wages, and computes what is left for the household's consumption, as
# main5's goods market clearing (section "Stationary equilibrium"):
#     C = Y - I - Psi
#   Y   = output of all producing firms
#   I   = delta * K: with a constant capital stock, investment replaces depreciation; the capital that liquidations
#         destroy is in Psi, as in main5
#   Psi = fixed operating costs + entry costs + profits lost in liquidation + capital lost in liquidation
#         + fixed (zeta) and variable ((1 - phi_c_hh) V) reorganization costs
# and from the household's side (main5's budget constraint): C = w N + T + net interest on the loans to firms.
# Labor supply is fixed, so the steady-state welfare gain in consumption units is the change in C.
#
# Timing. A firm that continues in state (x, e) with (k', b') produces next period with k' and the new productivity e',
# as in the law of motion (x' uses e'). sumSS and stat_dist measure output and labor with k' and the current e, so their
# firm mass clears a slightly different labor market. Here every aggregate is rescaled so that labor demand under the
# law-of-motion timing equals the labor supply of 10,000 ("labor scale" below); Y/L = w/nu under both timings.
#
# Entry costs. Free entry in the code is c_e = beta * E[V(0, e)] over all potential entrants (EntryValue), so every
# potential entrant pays c_e, also those who draw a low e and leave: entry costs are m * c_e.

const MODEL_DIR = dirname(@__DIR__)
const MODEL_FILE = joinpath(MODEL_DIR, get(ENV, "AGG_MODEL", "Baseline_noise60.jl"))
let code = read(MODEL_FILE, String)
    i = findfirst("########## " * "CALIBRATION", code)
    task_local_storage(() -> include_string(Main, code[1:first(i)-1], MODEL_FILE), :SOURCE_PATH, MODEL_FILE)
end
include(joinpath(MODEL_DIR, "Functions", "Results.jl"))

# welfare(E): the resource and income accounts of the equilibrium E (from solve_eq or ge_eq), per period.
# It follows each continuing firm through next period's events exactly as FirmOptim's Fmat does: production at e', then
# repayment, or default (chosen at the x' grid point, or the exogenous shock), then liquidation or reorganization.
function welfare(E)
    S, mu, p, w = E.SumPol, vec(E.mu), E.params, E.wage
    n = size(S, 1)                                  # the last row is the default state
    x_size, e_size = E.gridsize.x_size, E.gridsize.e_size
    alpha, nu, delta, pc = p.alpha, p.nu, p.delta, p.pc
    x_grid = S[1:x_size, 1]                         # x varies fastest in the state index
    e_vals = exp.(E.e_chain.state_values)
    P = E.e_chain.p
    V = S[:, 18]

    fn_L(k, e) = (nu * e * k^alpha / w)^(1 / (1 - nu))
    fn_Y(k, e) = e * k^alpha * fn_L(k, e)^nu
    # FirmOptim's liquidation probability in default
    fn_Gam(val, k, zeta_R) = p.sigma_v == 0 ? Float64(p.phi_a * (1 - delta) * k >= p.phi_c_hh * val - zeta_R) :
                             p.phi_c_hh * val <= 0 ? 1.0 : (p.phi_a * (1 - delta) * k + zeta_R) <= 0 ? 0.0 :
                             cdf(Normal(), log((p.phi_a * (1 - delta) * k + zeta_R) / (p.phi_c_hh * val)) / p.sigma_v)
    # V at cash on hand x and productivity index e_i, interpolated as in FirmOptim's q update (nearest edge off the grid)
    function V_at(x, e_i)
        base = (e_i - 1) * x_size
        j = searchsortedlast(x_grid, x)
        1 <= j < x_size || return V[(j <= 0 ? 1 : x_size) + base]
        t = (x - x_grid[j]) / (x_grid[j+1] - x_grid[j])
        (1 - t) * V[j + base] + t * V[j + 1 + base]
    end

    a = Dict{String, Float64}(k => 0.0 for k in
        ["firms", "K", "Y", "L", "exit mass", "exit value", "x - k'", "dividends", "new loans", "repaid",
         "default mass", "default mass, chosen", "default mass, borrowers", "liq mass, borrowers", "liq mass", "liq K", "liq profits", "reo mass", "reo V", "reo zeta",
         "grid cash - true cash", "mass on default states", "cash on default states"])
    arrivals = zeros(n)                             # next period's mass by state, to check against mu

    for s in 1:n-1
        m_s = mu[s]
        m_s == 0 && continue
        x = S[s, 1]
        if S[s, 7] != 0                             # a state that chooses default: stat_dist removes it (no mass lands here)
            a["mass on default states"] += m_s; a["cash on default states"] += m_s * x
            continue
        elseif S[s, 6] != 0                         # voluntary exit: the household receives x
            a["exit mass"] += m_s; a["exit value"] += m_s * x
            continue
        end
        k, b, q = S[s, 3], S[s, 4], S[s, 9]
        large = (p.size_by_debt == 1 ? b : k) >= p.Fcut
        zeta_R = large ? p.zeta_Rl : p.zeta_Rs
        pdef_exo = large ? p.pdef_exo_l : p.pdef_exo_s
        e_i = (s - 1) ÷ x_size + 1

        a["firms"] += m_s; a["K"] += m_s * k
        a["x - k'"] += m_s * (x - k); a["dividends"] += m_s * S[s, 13]; a["new loans"] += m_s * q * b

        for e_n in 1:e_size
            pt = P[e_i, e_n]
            pt == 0 && continue
            base = (e_n - 1) * x_size
            y = fn_Y(k, e_vals[e_n])
            profit = (1 - nu) * y - pc
            x_next = profit + (1 - delta) * k - b   # cash on hand after repaying b'
            a["Y"] += m_s * pt * y; a["L"] += m_s * pt * fn_L(k, e_vals[e_n])

            val_reo = V_at(x_next + b, e_n)         # reorganization: a fresh start at x' + b'
            liq = fn_Gam(val_reo, k, zeta_R)
            x_close, cw, x_far = landing(x_next, x_grid)
            for (x_i, wgt) in ((x_close, cw), (x_far, 1 - cw))
                (x_i == 0 || wgt == 0) && continue
                chosen = S[x_i + base, 7]
                p_def = chosen + (1 - chosen) * pdef_exo
                pm = m_s * pt * wgt

                # repays and stays at this grid point
                arrivals[x_i + base] += pm * (1 - p_def)
                a["repaid"] += pm * (1 - p_def) * b
                a["grid cash - true cash"] += pm * (1 - p_def) * (x_grid[x_i] - x_next)
                # defaults
                a["default mass"] += pm * p_def; a["default mass, chosen"] += pm * chosen
                a["liq mass"] += pm * p_def * liq
                if b != 0
                    a["default mass, borrowers"] += pm * p_def; a["liq mass, borrowers"] += pm * p_def * liq
                end
                a["liq K"] += pm * p_def * liq * k
                a["liq profits"] += pm * p_def * liq * profit
                reo = pm * p_def * (1 - liq)
                a["reo mass"] += reo; a["reo V"] += reo * val_reo; a["reo zeta"] += reo * zeta_R
                if reo > 0                          # reorganized: lands around x' + b'
                    r_close, rw, r_far = landing(x_next + b, x_grid)
                    arrivals[r_close + base] += reo * rw
                    a["grid cash - true cash"] += reo * rw * (x_grid[r_close] - (x_next + b))
                    if r_far != 0
                        arrivals[r_far + base] += reo * (1 - rw)
                        a["grid cash - true cash"] += reo * (1 - rw) * (x_grid[r_far] - (x_next + b))
                    end
                end
            end
        end
    end

    # check: the events above reproduce the stationary distribution, mu = entrants that stay + arrivals
    entrants = E.m .* (1 .- vec(E.xpol)) .* vec(E.f0)
    law_gap = maximum(abs.(arrivals[1:n-1] .+ entrants[1:n-1] .- mu[1:n-1])) / sum(mu[1:n-1])

    sc = 10000 / a["L"]                             # labor scale: labor demand at the law-of-motion timing = 10,000
    r = Dict{String, Float64}()
    r["Y"] = sc * a["Y"]
    r["I = delta K"] = sc * delta * a["K"]
    r["Psi: fixed operating costs"] = sc * pc * a["firms"]
    r["Psi: entry costs"] = sc * E.m * E.c_e
    r["Psi: capital lost in liquidation"] = sc * (1 - p.phi_a) * (1 - delta) * a["liq K"]
    r["Psi: profits lost in liquidation"] = sc * a["liq profits"]
    r["Psi: reorganization, fixed (zeta)"] = sc * a["reo zeta"]
    r["Psi: reorganization, variable"] = sc * (1 - p.phi_c_hh) * a["reo V"]
    r["Psi"] = sum(v for (k, v) in r if startswith(k, "Psi: "))
    r["C = Y - I - Psi"] = r["Y"] - r["I = delta K"] - r["Psi"]
    r["C if reorganization costs are transfers"] = r["C = Y - I - Psi"] + r["Psi: reorganization, fixed (zeta)"] + r["Psi: reorganization, variable"]

    # the household's side: wage income, and everything it receives from and pays to firms and the lender.
    # The lender's recoveries in default cancel out: what the household pays in a reorganization (or leaves to the
    # lender in a liquidation) returns to it as the lender's owner.
    r["HH: wage income w N"] = sc * nu * a["Y"]
    r["HH: dividends"] = sc * a["dividends"]
    r["HH: exit values"] = sc * a["exit value"]
    r["HH: liquidation proceeds"] = sc * p.phi_a * (1 - delta) * a["liq K"]
    r["HH: loans repaid - new loans"] = sc * (a["repaid"] - a["new loans"])
    r["HH: entry costs"] = -r["Psi: entry costs"]
    r["HH: reorganization costs"] = -(r["Psi: reorganization, fixed (zeta)"] + r["Psi: reorganization, variable"])
    r["HH: total"] = sum(v for (k, v) in r if startswith(k, "HH: "))
    # The x-grid lottery moves cash: a firm lands on a grid point above or below its true x', and the two points can
    # differ in default. "grid cash" is that created cash; the household account holds it, the resource account does not.
    r["check: grid cash created"] = sc * (a["grid cash - true cash"] - a["cash on default states"])
    r["check: HH total - grid cash - C"] = r["HH: total"] - r["check: grid cash created"] - r["C = Y - I - Psi"]
    r["check: law of motion gap (share of mass)"] = law_gap
    r["check: mass on default states"] = sc * a["mass on default states"]

    r["labor scale"] = sc
    r["producing firms"] = sc * a["firms"]
    r["potential entrants"] = sc * E.m
    r["K"] = sc * a["K"]
    r["default rate, %"] = 100 * a["default mass"] / a["firms"]
    r["liquidations per 100 firms"] = 100 * a["liq mass"] / a["firms"]
    r["reorganizations per 100 firms"] = 100 * a["reo mass"] / a["firms"]
    r["voluntary exits per 100 firms"] = 100 * a["exit mass"] / a["firms"]
    # liquidations as a share of defaults (a flow; sumSS's avg_liqprob is the firm average of the liquidation probability)
    r["liquidated share of defaults, %"] = 100 * a["liq mass"] / a["default mass"]
    r["liquidated share of defaults, borrowers, %"] = 100 * a["liq mass, borrowers"] / a["default mass, borrowers"]
    r["borrowers' share of defaults, %"] = 100 * a["default mass, borrowers"] / a["default mass"]
    r["chosen share of defaults, %"] = 100 * a["default mass, chosen"] / a["default mass"]
    r["avg_liqprob (sumSS), %"] = 100 * E.ss["avg_liqprob", "values"]
    r["wage"] = w
    return r
end

reset_params!()
@time BaseEq = solve_eq()
@time Reform2 = ge_eq((zeta_Rl = parameters().zeta_Rl / 2,), BaseEq.c_e)
@time Reform3 = ge_eq((zeta_Rl = parameters().zeta_Rl / 3,), BaseEq.c_e)

wb, w2, w3 = welfare(BaseEq), welfare(Reform2), welfare(Reform3)
C0 = wb["C = Y - I - Psi"]
num(v; d = 2) = lpad(round(v, digits = d), 13)

# Levels, and each item's change as a percent of baseline consumption: the rows of C = Y - I - Psi add up to the change in C
println("\nResources, per period (labor supply 10,000). 'pp of C': change in percent of baseline C; minus the costs' changes add up to C's")
println(rpad("", 42), lpad("baseline", 13), lpad("zeta/2", 13), lpad("zeta/3", 13), lpad("z/2 pp of C", 13), lpad("z/3 pp of C", 13))
for k in ["Y", "I = delta K", "Psi: fixed operating costs", "Psi: entry costs", "Psi: capital lost in liquidation",
          "Psi: profits lost in liquidation", "Psi: reorganization, fixed (zeta)", "Psi: reorganization, variable", "Psi",
          "C = Y - I - Psi", "C if reorganization costs are transfers"]
    println(rpad(k, 42), num(wb[k]), num(w2[k]), num(w3[k]), num(100 * (w2[k] - wb[k]) / C0; d = 3), num(100 * (w3[k] - wb[k]) / C0; d = 3))
end

println("\nHousehold's side (main5's budget: C = w N + T + net interest)")
println(rpad("", 42), lpad("baseline", 13), lpad("zeta/2", 13), lpad("zeta/3", 13), lpad("z/2 pp of C", 13), lpad("z/3 pp of C", 13))
for k in ["HH: wage income w N", "HH: dividends", "HH: exit values", "HH: liquidation proceeds", "HH: loans repaid - new loans",
          "HH: entry costs", "HH: reorganization costs", "HH: total", "check: grid cash created", "check: HH total - grid cash - C"]
    println(rpad(k, 42), num(wb[k]), num(w2[k]), num(w3[k]), num(100 * (w2[k] - wb[k]) / C0; d = 3), num(100 * (w3[k] - wb[k]) / C0; d = 3))
end

println("\nOther")
println(rpad("", 42), lpad("baseline", 13), lpad("zeta/2", 13), lpad("zeta/3", 13), lpad("z/2 %", 13), lpad("z/3 %", 13))
for k in ["wage", "labor scale", "producing firms", "potential entrants", "K", "default rate, %", "liquidations per 100 firms",
          "reorganizations per 100 firms", "voluntary exits per 100 firms", "liquidated share of defaults, %",
          "liquidated share of defaults, borrowers, %", "borrowers' share of defaults, %", "chosen share of defaults, %",
          "avg_liqprob (sumSS), %", "check: mass on default states"]
    pc_(x) = wb[k] == 0 ? NaN : 100 * (x / wb[k] - 1)
    println(rpad(k, 42), num(wb[k]; d = 5), num(w2[k]; d = 5), num(w3[k]; d = 5), num(pc_(w2[k]); d = 2), num(pc_(w3[k]); d = 2))
end

println("check: law of motion gap (largest state gap, share of mass): ", wb["check: law of motion gap (share of mass)"], ", ",
        w2["check: law of motion gap (share of mass)"], ", ", w3["check: law of motion gap (share of mass)"])
println("\nWelfare gain (change in steady-state consumption): zeta/2 ", round(100 * (w2["C = Y - I - Psi"] / C0 - 1), digits = 3),
        "%, zeta/3 ", round(100 * (w3["C = Y - I - Psi"] / C0 - 1), digits = 3), "%")
println("Y/L gain (sumSS):                                  zeta/2 ", round(100 * (Reform2.ss["YtoL", "values"] / BaseEq.ss["YtoL", "values"] - 1), digits = 3),
        "%, zeta/3 ", round(100 * (Reform3.ss["YtoL", "values"] / BaseEq.ss["YtoL", "values"] - 1), digits = 3), "%")

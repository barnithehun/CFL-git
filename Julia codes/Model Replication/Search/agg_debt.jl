############ Aggregate vs firm-average debt ratios: baseline and zeta/2 reform #########
# Loads a model file up to its CALIBRATION block (default Baseline_noise60.jl, ENV AGG_MODEL to change), solves the
# baseline and the zeta/2 reform at its GE wage, and compares debt to collateral and debt to output:
#   firm average   = mean over firms of b'/k' (sumSS's measure: q >= 0.5, b'/k' <= d2c_cut) or b'/y
#   economy wide   = sum of b' over sum of k' (or y), all continuing firms, and the same among borrowers only
# Small/large split as sumSSsme (borrowers, b' < Fcut is small). b' is the face value of debt, as in sumSS.

const MODEL_DIR = dirname(@__DIR__)
const MODEL_FILE = joinpath(MODEL_DIR, get(ENV, "AGG_MODEL", "Baseline_noise60.jl"))
let code = read(MODEL_FILE, String)
    i = findfirst("########## " * "CALIBRATION", code)
    task_local_storage(() -> include_string(Main, code[1:first(i)-1], MODEL_FILE), :SOURCE_PATH, MODEL_FILE)
end
include(joinpath(MODEL_DIR, "Functions", "Results.jl"))

function debt_ratios(eq)
    S, mu = eq.SumPol, vec(eq.mu)
    n = size(S, 1)
    Fcut, cut = eq.params.Fcut, eq.params.d2c_cut
    out = Dict{String, Float64}()
    groups = ["All firms" => (i -> true), "Small" => (i -> S[i, 4] != 0 && S[i, 4] < Fcut),
              "Large" => (i -> S[i, 4] != 0 && S[i, 4] >= Fcut)]
    for (g, ing) in groups
        K = B = Y = Kb = Yb = 0.0                              # Kb, Yb: borrowers only
        fa = fw = fy = fyw = 0.0                               # firm averages
        for i in 1:n-1
            ing(i) || continue
            k, b, y, m = S[i, 3], S[i, 4], S[i, 11], mu[i]
            K += m * k; B += m * b; Y += m * y
            if b != 0
                Kb += m * k; Yb += m * y
            end
            if k != 0 && S[i, 9] >= 0.5 && b / k <= cut        # sumSS's b'/k' measure
                fa += m * b / k; fw += m
            end
            if y != 0 && b != 0                                # borrowers' b'/y
                fy += m * b / y; fyw += m
            end
        end
        out[g * ": B/K firm avg"] = fa / fw
        out[g * ": B/K aggregate, all"] = B / K
        out[g * ": B/K aggregate, borrowers"] = B / Kb
        out[g * ": B/Y firm avg (borrowers)"] = fy / fyw
        out[g * ": B/Y aggregate, all"] = B / Y
        out[g * ": B/Y aggregate, borrowers"] = B / Yb
        out[g * ": total B"] = B
        out[g * ": total K"] = K
        out[g * ": total Y"] = Y
    end
    out["borrower share"] = sum(mu[i] for i in 1:n-1 if S[i, 4] != 0) / sum(mu)
    out["mass"] = sum(mu)

    # Allocation of debt and capital across borrowers (k' > 0, b' != 0). e = productivity level (SumPol col 2),
    # y/k = average product of capital (MPK = alpha y/k, same alpha for all firms)
    bor = [i for i in 1:n-1 if S[i, 4] != 0 && S[i, 3] > 0 && mu[i] > 0]
    m, e, k, b, y = mu[bor], S[bor, 2], S[bor, 3], S[bor, 4], S[bor, 11]
    wmean(x, w) = sum(w .* x) / sum(w)
    wsd(x, w) = sqrt(wmean((x .- wmean(x, w)) .^ 2, w))
    le = log.(e)
    out["A1 mean e, firm-weighted"] = wmean(e, m)
    out["A2 mean e, debt-weighted"] = wmean(e, m .* b)
    out["A3 mean e, capital-weighted"] = wmean(e, m .* k)
    out["A4 OP cov debt share x e (A2-A1)"] = wmean(e, m .* b) - wmean(e, m)
    out["A5 OP cov debt share x log e"] = wmean(le, m .* b) - wmean(le, m)
    out["A6 OP cov capital share x log e"] = wmean(le, m .* k) - wmean(le, m)
    out["A7 sd log(y/k), firm-weighted"] = wsd(log.(y ./ k), m)
    out["A8 sd log(y/k), capital-weighted"] = wsd(log.(y ./ k), m .* k)
    out["A9 small firms' share of debt"] = sum(m[j] * b[j] for j in eachindex(b) if b[j] < Fcut) / sum(m .* b)
    return out
end

reset_params!()
@time BaseEq = solve_eq()
@time Reform2 = ge_eq((zeta_Rl = parameters().zeta_Rl / 2,), BaseEq.c_e)
@time Reform3 = ge_eq((zeta_Rl = parameters().zeta_Rl / 3,), BaseEq.c_e)
for (nm, R) in ["zeta/2" => Reform2, "zeta/3" => Reform3]
    println(nm, ": wage ", R.wage, "; sumSS avg_b2a ", BaseEq.ss["avg_b2a", "values"], " -> ", R.ss["avg_b2a", "values"],
            "; Y/L ", round(100 * (R.ss["YtoL", "values"] / BaseEq.ss["YtoL", "values"] - 1), digits = 2), "%")
end

rb, r2, r3 = debt_ratios(BaseEq), debt_ratios(Reform2), debt_ratios(Reform3)
pct(a, b) = round(100 * (a / b - 1), digits = 2)
println(rpad("measure", 44), lpad("baseline", 12), lpad("zeta/2", 12), lpad("zeta/3", 12), lpad("chg2 %", 9), lpad("chg3 %", 9))
for k in sort(collect(keys(rb)))
    println(rpad(k, 44), lpad(round(rb[k], digits = 4), 12), lpad(round(r2[k], digits = 4), 12), lpad(round(r3[k], digits = 4), 12),
            lpad(pct(r2[k], rb[k]), 9), lpad(pct(r3[k], rb[k]), 9))
end

# liquidation probability (conditional on default): all firms (sumSS) and small/large (sumSSsme)
println(rpad("liquidation probability", 44), lpad("baseline", 12), lpad("zeta/2", 12), lpad("zeta/3", 12))
for (lab, f) in ["All firms" => (E -> E.ss["avg_liqprob", "values"]),
                 "Small" => (E -> E.sme["liqprob", "SME Values"]), "Large" => (E -> E.sme["liqprob", "LE Values"])]
    println(rpad(lab, 44), lpad(round(f(BaseEq), digits = 4), 12), lpad(round(f(Reform2), digits = 4), 12), lpad(round(f(Reform3), digits = 4), 12))
end

# cash-flow debt: CF share = CF debt / total debt (debt-weighted); CF reliance = firm average of tau (borrowers)
println(rpad("cash-flow debt", 44), lpad("baseline", 12), lpad("zeta/2", 12), lpad("zeta/3", 12))
for (lab, f) in ["CF share, all" => (E -> E.ss["CF share", "values"]), "CF reliance, all" => (E -> E.ss["CF reliance", "values"]),
                 "CF share, small" => (E -> E.sme["CF share", "SME Values"]), "CF share, large" => (E -> E.sme["CF share", "LE Values"]),
                 "CF reliance, small" => (E -> E.sme["CF reliance", "SME Values"]), "CF reliance, large" => (E -> E.sme["CF reliance", "LE Values"])]
    println(rpad(lab, 44), lpad(round(f(BaseEq), digits = 4), 12), lpad(round(f(Reform2), digits = 4), 12), lpad(round(f(Reform3), digits = 4), 12))
end

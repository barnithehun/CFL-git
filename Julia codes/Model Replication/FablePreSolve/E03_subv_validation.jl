############ E03: the targeted reform against the Subchapter V evidence ############
# Hotchkiss, Iverson and Zheng (2024): under Subchapter V eligible firms confirm a plan 21 pp more often (base 33%),
# filings do not rise, creditor recoveries and post-bankruptcy survival are not lower. Here: the targeted reform of
# the RESULTS section (small firms' zeta cut so that their liquidation probability falls by 13 pp at the baseline wage),
# in GE, and the universal zeta/3 for comparison. Small firms: borrowers with b' < Fcut, the model's size split.
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("E03")

B = base_eq()
SBRA = cached(() -> (s = find_zeta_Rs(B; drop = 0.13); (zeta_Rs = s.zeta_Rs, liqprob = s.liqprob, target = s.target)), "sbra_zeta_Rs")
Tg = cached(() -> ge((zeta_Rs = SBRA.zeta_Rs,)), "targeted")
R3 = reform_eq(3)
println("targeted: zeta_Rs = ", SBRA.zeta_Rs, " (", SBRA.zeta_Rs / B.params.zeta_Rl, " of zeta), GE wage ", Tg.wage)

function group_stats(E, small::Bool)
    ev = events(E); sv = statevecs(E, ev); K = kernel(ev); c = cont_mass(E); S = E.SumPol; p = E.params
    b_s, k_s = S[1:ev.n, 4], S[1:ev.n, 3]
    insize = ((p.size_by_debt == 1 ? b_s : k_s) .< p.Fcut) .== small
    g = (b_s .> 0) .& insize                               # the group's states: borrowers of this size
    w = c[ev.s] .* ev.prob
    rows = g[ev.s]
    rec = recoveries(E, ev)
    d, l, r = w .* ev.pdef, w .* ev.pdef .* ev.liq, w .* ev.pdef .* (1 .- ev.liq)
    # reorganized firms of the group, followed for 3 years, against the group's average firm
    d_reo = K.Treo * (c .* g)
    f_reo = follow(E, K, sv, d_reo; T = 3)
    f_avg = follow(E, K, sv, K.T * (c .* g); T = 3)
    st = state_stats(E, c .* g)
    return ["default rate, % of the group's borrowers" => 100 * sum(d[rows]) / sum(c[g]),
            "reorganized share of defaults, %" => 100 * sum(r[rows]) / sum(d[rows]),
            "liquidation prob., firm average, %" => 100 * E.sme["liqprob", small ? "SME Values" : "LE Values"],
            "lender recovery in default, % of face" => 100 * (sum(l[rows] .* rec.liq[rows]) + sum(r[rows] .* rec.reo[rows])) / sum(d[rows] .* ev.b[rows]),
            "lender recovery in reorganization, %" => 100 * sum(r[rows] .* rec.reo[rows]) / sum(r[rows] .* ev.b[rows]),
            "interest rate, %" => st.rate,
            "CF reliance" => st.tau,
            "reorganized firms: operating after 3 years, %" => 100 * f_reo.alive[end],
            "reorganized firms: default again within 3 years, %" => 100 * f_reo.cdef[end],
            "average firm: operating after 3 years, %" => 100 * f_avg.alive[end],
            "average firm: default within 3 years, %" => 100 * f_avg.cdef[end]]
end

for (small, file) in [(true, "small_firms.csv"), (false, "large_firms.csv")]
    rows = [group_stats(E, small) for E in (B, Tg, R3)]
    save_table(dir, file, ["baseline", "targeted", "zeta/3 all"], [first(rows[1][i]) => [last(r[i]) for r in rows] for i in eachindex(rows[1])])
end
ob = outcomes(B)
save_table(dir, "aggregates.csv", SWEEP_COLS, ["targeted" => sweep_row(outcomes(Tg), ob), "zeta/3 all" => sweep_row(outcomes(R3), ob)])

s = [group_stats(E, true) for E in (B, Tg)]
reo(i) = last(s[i][2])
plt = small_plot(title = "Reorganized share of small firms' defaults, %", ylabel = "%", ylims = (0, 80))
bar!(plt, [1, 2], [reo(1), reo(2)], bar_width = 0.6, label = "model", xticks = ([1, 2, 4, 5], ["baseline", "targeted", "data: Ch. 11", "data: Subch. V"]))
bar!(plt, [4, 5], [32.9, 32.9 + 21], bar_width = 0.6, label = "Hotchkiss et al. (plan confirmed)")
savefig(plt, joinpath(dir, "reorganized_share.png"))

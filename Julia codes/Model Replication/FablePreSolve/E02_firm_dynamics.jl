############ E02: firm dynamics and lender recoveries ############
# Entry, exit and default rates, how defaults are resolved, and what lenders recover per unit of face value in
# liquidation and in reorganization, for mostly asset-based (tau < 0.5) and mostly CF-based (tau >= 0.5) borrowers.
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("E02")

function dynamics(E)
    ev = events(E); sv = statevecs(E, ev); c = cont_mass(E); S = E.SumPol
    w = c[ev.s] .* ev.prob                       # mass of each event row
    rec = recoveries(E, ev)
    d, l, r = w .* ev.pdef, w .* ev.pdef .* ev.liq, w .* ev.pdef .* (1 .- ev.liq)
    bor, cf = ev.b .> 0, ev.tau .>= 0.5
    firms = sum(c)
    borrowers = sum(c[S[1:ev.n, 4] .> 0])
    exits = sv.exit' * vec(E.mu)[1:end-1]
    rr(mask, flow, got) = 100 * sum(flow[mask] .* got[mask]) / sum(flow[mask] .* ev.b[mask])
    return ["entry rate, % of firms" => 100 * E.m * sum(entrant_dist(E)) / firms,
            "potential entrants that stay, %" => 100 * sum(entrant_dist(E)),
            "exit rate, % (voluntary + liquidated)" => 100 * (exits + sum(l)) / firms,
            "  voluntary exits, %" => 100 * exits / firms,
            "  liquidations, %" => 100 * sum(l) / firms,
            "default rate, % of firms" => 100 * sum(d) / firms,
            "default rate, % of borrowers" => 100 * sum(d[bor]) / borrowers,
            "defaults without debt, % of defaults" => 100 * sum(d[.!bor]) / sum(d),
            "liquidated share of defaults, borrowers, %" => 100 * sum(l[bor]) / sum(d[bor]),
            "liquidation prob., firm average (sumSS), %" => 100 * E.ss["avg_liqprob", "values"],
            "recovery in liquidation, % of face value" => rr(bor, l, rec.liq),
            "recovery in reorganization, %" => rr(bor, r, rec.reo),
            "recovery in default, %" => 100 * (sum(l[bor] .* rec.liq[bor]) + sum(r[bor] .* rec.reo[bor])) / sum(d[bor] .* ev.b[bor]),
            "AB borrowers: recovery in liquidation, %" => rr(bor .& .!cf, l, rec.liq),
            "AB borrowers: recovery in reorganization, %" => rr(bor .& .!cf, r, rec.reo),
            "CF borrowers: recovery in liquidation, %" => rr(bor .& cf, l, rec.liq),
            "CF borrowers: recovery in reorganization, %" => rr(bor .& cf, r, rec.reo),
            "CF borrowers' share of defaulted debt, %" => 100 * sum(d[bor .& cf] .* ev.b[bor .& cf]) / sum(d[bor] .* ev.b[bor]),
            "interest rate, firm average, %" => E.ss["avg intrate", "values"]]
end

rows = [dynamics(base_eq()), dynamics(reform_eq(2)), dynamics(reform_eq(3))]
tab = [first(rows[1][i]) => [last(r[i]) for r in rows] for i in eachindex(rows[1])]
save_table(dir, "firm_dynamics.csv", ["baseline", "zeta/2", "zeta/3"], tab)

val(label, j) = last(tab[findfirst(t -> first(t) == label, tab)])[j]
labels = ["AB, liquidation", "AB, reorganization", "CF, liquidation", "CF, reorganization"]
keys_ = ["AB borrowers: recovery in liquidation, %", "AB borrowers: recovery in reorganization, %",
         "CF borrowers: recovery in liquidation, %", "CF borrowers: recovery in reorganization, %"]
plt = small_plot(title = "Lender recovery, % of face value", xrotation = 20, bottom_margin = 7Measures.mm, ylabel = "%")
bar!(plt, (1:4) .- 0.2, [val(k, 1) for k in keys_], bar_width = 0.38, label = "baseline", xticks = (1:4, labels))
bar!(plt, (1:4) .+ 0.2, [val(k, 3) for k in keys_], bar_width = 0.38, label = "zeta/3")
savefig(plt, joinpath(dir, "recoveries.png"))

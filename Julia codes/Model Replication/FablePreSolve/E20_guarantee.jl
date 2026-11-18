############ E20: a reorganization subsidy against a credit guarantee, per unit of public money ############
# Subsidy: the state pays the share s of zeta. Firms, lenders and the court face (1 - s) zeta, the economy still uses
# up zeta per reorganization. Fiscal cost: s * zeta per reorganization.
# Guarantee: the state tops up what a lender recovers on CF debt in a liquidation, from the share kappa of the
# liquidation value to kappa' (kappa' = 1: CF debt recovers as much as asset-based debt). Fiscal cost: the top-up.
# Both are paid with lump-sum taxes, so the fiscal cost is a transfer; consumption counts real resources only.
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("E20")
B = base_eq(); ob = outcomes(B); zeta, kappa0 = B.params.zeta_Rl, B.params.kappa
CB = ob["C"]

# the guarantee's cost: liquidations' recoveries at kappa' minus at the baseline kappa, per period, labor 10,000
function guarantee_cost(E)
    ev = events(E); c = cont_mass(E); p = E.params
    sc = 10000 / (statevecs(E, ev).L' * c)
    l = c[ev.s] .* ev.prob .* ev.pdef .* ev.liq
    Pl = max.(0, p.phi_af .* (1 - p.delta) .* ev.k .- p.zeta_L)
    pay(kap) = min.(ev.b, (1 .- ev.tau) .* Pl .+ ev.tau .* kap .* Pl)
    return sc * sum(l .* (pay(p.kappa) .- pay(kappa0)))
end

subs = [0.25, 0.5, 2/3, 0.9]
kaps = [0.3, 0.5, 0.75, 1.0]
res = cached("subsidy_guarantee") do
    sub = map(subs) do s
        o = outcomes(ge((zeta_Rl = (1 - s) * zeta,)))
        cost = s / (1 - s) * o["reo zeta"]                  # o["reo zeta"] is (1 - s) zeta times the reorganizations
        (o = o, cost = cost, C = o["C"] - cost)             # the economy still pays the subsidized part
    end
    gua = map(kaps) do kap
        E = ge((kappa = kap,))
        o = outcomes(E)
        (o = o, cost = guarantee_cost(E), C = o["C"])
    end
    (sub = sub, gua = gua)
end
cols = ["Y/L %", "C %", "fiscal cost, % of C", "C gain per unit cost", "CF reliance", "liq share", "rate small", "default rate"]
row(r) = [pct(r.o["wage"], ob["wage"]), pct(r.C, CB), 100 * r.cost / CB, (r.C - CB) / r.cost, r.o["CF reliance"], 100 * r.o["liq share"],
          r.o["rate small"], r.o["default rate"]]
save_table(dir, "subsidy_vs_guarantee.csv", cols,
           [["baseline" => [0.0, 0.0, 0.0, NaN, ob["CF reliance"], 100 * ob["liq share"], ob["rate small"], ob["default rate"]]];
            ["subsidy $(round(Int, 100s))% of zeta" => row(r) for (s, r) in zip(subs, res.sub)];
            ["guarantee kappa $k" => row(r) for (k, r) in zip(kaps, res.gua)]])

plt = small_plot(title = "Consumption gain against fiscal cost, % of baseline C", xlabel = "fiscal cost", ylabel = "consumption gain")
plot!(plt, [100 * r.cost / CB for r in res.sub], [pct(r.C, CB) for r in res.sub], marker = :circle, label = "reorganization subsidy")
plot!(plt, [100 * r.cost / CB for r in res.gua], [pct(r.C, CB) for r in res.gua], marker = :square, label = "credit guarantee")
hline!(plt, [0], color = :black, label = "")
savefig(plt, joinpath(dir, "subsidy_vs_guarantee.png"))

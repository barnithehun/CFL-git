############ E12: is the court's liquidation rule the best one? ############
# The court liquidates when the liquidation value beats the reorganization value net of costs: right after the
# default, but it ignores that a higher chance of liquidation makes CF debt dearer before the default. liq_wedge tilts
# the rule: the court acts as if the reorganization value were e^wedge times higher (wedge > 0: fewer liquidations).
# Each wedge is solved in GE; the best rule is the one with the highest steady-state consumption.
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("E12")
B = base_eq(); ob = outcomes(B); zeta = B.params.zeta_Rl

wedges = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
res = cached("wedge_sweep") do
    [w == 0.0 ? ob : outcomes(ge((liq_wedge = w,); bracket = (0.8, 1.3))) for w in wedges]
end
wedges3 = [-1.0, 0.0, 1.0, 2.0]
res3 = cached("wedge_sweep_zeta3") do
    [w == 0.0 ? outcomes(reform_eq(3)) : outcomes(ge((liq_wedge = w, zeta_Rl = zeta / 3); bracket = (0.8, 1.3))) for w in wedges3]
end
save_table(dir, "wedge.csv", SWEEP_COLS,
           [["wedge $w" => sweep_row(o, ob) for (w, o) in zip(wedges, res)]; ["zeta/3, wedge $w" => sweep_row(o, ob) for (w, o) in zip(wedges3, res3)]])
i = argmax([o["C"] for o in res])
println("\nbest wedge for consumption on the grid: ", wedges[i], " (C ", round(pct(res[i]["C"], ob["C"]), digits = 3), "%, Y/L ",
        round(pct(res[i]["wage"], ob["wage"]), digits = 3), "%); main5 entry costs: ", wedges[argmax([o["C main5"] for o in res])])

p1 = small_plot(title = "Change from the baseline court, %", xlabel = "wedge toward reorganization (log points)")
plot!(p1, wedges, [pct(o["wage"], ob["wage"]) for o in res], marker = :circle, label = "Y/L")
plot!(p1, wedges, [pct(o["C"], ob["C"]) for o in res], marker = :square, label = "consumption")
plot!(p1, wedges, [pct(o["C main5"], ob["C main5"]) for o in res], marker = :diamond, linestyle = :dash, label = "consumption, main5 entry costs")
hline!(p1, [0], color = :black, label = "")
p2 = small_plot(title = "Defaults and debt", xlabel = "wedge toward reorganization (log points)")
plot!(p2, wedges, [o["liq share"] for o in res], marker = :circle, label = "liquidated share of defaults")
plot!(p2, wedges, [o["CF reliance"] for o in res], marker = :square, label = "CF reliance")
savefig(plot(p1, p2, layout = (1, 2), size = (900, 340), dpi = 110, bottom_margin = 5Measures.mm), joinpath(dir, "wedge.png"))

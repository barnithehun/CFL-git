############ E12b: the court's liquidation rule without valuation noise ############
# E12 on noise46 mixes two things: a court that ignores the effect of liquidations on loan prices, and a court that
# makes valuation errors, which a wedge toward reorganization partly undoes. Here the same sweep is run on the
# noiseless model: Baseline.jl's calibration (DRS 0.85, sigma_v = 0, main5's tie rule), set as overrides on Model_ext.jl,
# which is the same model code at these values. The court then liquidates exactly when L >= phi_c_hh V e^wedge - zeta.
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("E12b")

# parameters() of Baseline.jl (2026-11-06 full calibration, error 0.057, c_e 340.776)
N0 = (pc = 30.221705598739636, pdef_exo_l = 0.01623279813010318, phi_a = 0.4, phi_af = 0.4, phi_c_hh = 0.7176063676733062,
      Fcut = 533.6699231206312, phi_c = 0.21842637289640024, zeta_Rl = 1154.3576752772226, sigma_v = 0.0, tau_tie = 0.0)
reset_params!()
B0 = solve_eq(N0)
println("noiseless baseline: c_e ", B0.c_e, " (Baseline.jl: 340.776)")
ob = with_params(() -> outcomes(B0), N0)

wedges = [-1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
res = cached("wedge_sweep_noiseless") do
    map(wedges) do w
        w == 0.0 && return ob
        vals = merge(N0, (liq_wedge = w,))
        with_params(() -> outcomes(ge_eq(vals, B0.c_e; bracket = (0.8, 1.3))), vals)
    end
end
ref = cached("noiseless_zeta_reforms") do
    map((2, 3)) do dv
        vals = merge(N0, (zeta_Rl = N0.zeta_Rl / dv,))
        with_params(() -> outcomes(ge_eq(vals, B0.c_e; bracket = (0.8, 1.3))), vals)
    end
end
save_table(dir, "wedge_noiseless.csv", SWEEP_COLS,
           [["wedge $w" => sweep_row(o, ob) for (w, o) in zip(wedges, res)]; ["zeta/2, no wedge" => sweep_row(ref[1], ob), "zeta/3, no wedge" => sweep_row(ref[2], ob)]])
i = argmax([o["C"] for o in res])
println("\nbest wedge for consumption on the grid: ", wedges[i], " (C ", round(pct(res[i]["C"], ob["C"]), digits = 3), "%, Y/L ",
        round(pct(res[i]["wage"], ob["wage"]), digits = 3), "%); main5 entry costs: ", wedges[argmax([o["C main5"] for o in res])])

p1 = small_plot(title = "Noiseless court: change from wedge 0, %", xlabel = "wedge toward reorganization (log points)")
plot!(p1, wedges, [pct(o["wage"], ob["wage"]) for o in res], marker = :circle, label = "Y/L")
plot!(p1, wedges, [pct(o["C"], ob["C"]) for o in res], marker = :square, label = "consumption")
plot!(p1, wedges, [pct(o["C main5"], ob["C main5"]) for o in res], marker = :diamond, linestyle = :dash, label = "consumption, main5 entry costs")
hline!(p1, [0], color = :black, label = "")
p2 = small_plot(title = "Defaults and debt", xlabel = "wedge toward reorganization (log points)")
plot!(p2, wedges, [o["liq share"] for o in res], marker = :circle, label = "liquidated share of defaults")
plot!(p2, wedges, [o["CF reliance"] for o in res], marker = :square, label = "CF reliance")
savefig(plot(p1, p2, layout = (1, 2), size = (900, 340), dpi = 110, bottom_margin = 5Measures.mm), joinpath(dir, "wedge_noiseless.png"))

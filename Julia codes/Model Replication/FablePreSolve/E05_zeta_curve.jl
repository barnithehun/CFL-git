############ E05: the whole zeta curve ############
# The fixed reorganization cost from zero to prohibitive (liquidation only), each at its GE wage, and the perfect-credit
# benchmark: how much of the cost of credit frictions do bankruptcy costs account for?
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("E05")
B = base_eq(); ob = outcomes(B); zeta = B.params.zeta_Rl

mults = [0.0, 0.1, 0.25, 1/3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 20.0]
res = cached("zeta_curve") do
    map(mults) do m
        E = m == 1.0 ? B : ge((zeta_Rl = m * zeta,); bracket = (0.8, 1.3))
        outcomes(E)
    end
end
rows = [(m == 20.0 ? "zeta x 20 (liquidation only)" : "zeta x $(round(m, digits = 3))") => sweep_row(o, ob) for (m, o) in zip(mults, res)]

# perfect credit: negative dividends are free (perf_cred = 1), as in the RESULTS section
perf = cached("perfect_credit") do
    try
        outcomes(ge((perf_cred = 1.0,); bracket = (1.0, 1.6)))
    catch err
        println("perfect credit failed: ", err); nothing
    end
end
perf === nothing || push!(rows, "perfect credit" => sweep_row(perf, ob))
save_table(dir, "zeta_curve.csv", SWEEP_COLS, rows)

if perf !== nothing
    gap = pct(perf["wage"], ob["wage"])
    lo = pct(res[end]["wage"], ob["wage"])
    println("\nY/L: liquidation only ", round(lo, digits = 2), "%, zeta = 0 ", round(pct(res[1]["wage"], ob["wage"]), digits = 2),
            "%, perfect credit ", round(gap, digits = 2), "%")
    println("share of the baseline-to-perfect-credit gap closed: zeta/3 ", round(100 * pct(res[4]["wage"], ob["wage"]) / gap, digits = 1),
            "%, zeta = 0 ", round(100 * pct(res[1]["wage"], ob["wage"]) / gap, digits = 1), "%")
end

x = 1:length(mults)
labs = ["0", "0.1", "1/4", "1/3", "1/2", "3/4", "1", "1.5", "2", "3", "5", "20"]
p1 = small_plot(title = "Change from the baseline, %", xlabel = "zeta, multiple of the baseline", xticks = (x, labs))
plot!(p1, x, [pct(o["wage"], ob["wage"]) for o in res], marker = :circle, label = "Y/L")
plot!(p1, x, [pct(o["C"], ob["C"]) for o in res], marker = :square, label = "consumption")
hline!(p1, [0], color = :black, label = "")
p2 = small_plot(title = "Defaults and debt", xlabel = "zeta, multiple of the baseline", xticks = (x, labs))
plot!(p2, x, [o["liq share"] for o in res], marker = :circle, label = "liquidated share of defaults")
plot!(p2, x, [o["CF reliance"] for o in res], marker = :square, label = "CF reliance")
savefig(plot(p1, p2, layout = (1, 2), size = (900, 340), dpi = 110, bottom_margin = 5Measures.mm), joinpath(dir, "zeta_curve.png"))

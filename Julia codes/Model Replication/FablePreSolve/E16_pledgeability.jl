############ E16: does the reform matter more when assets are worth less in liquidation? ############
# phi_a is the resale value of capital in liquidation (Kermani and Ma: lower in industries with specific assets; a
# lower phi_a also stands for more intangible capital). For each phi_a the economy is solved in GE with the baseline
# entry cost and no refit, with zeta and with zeta/3. The reform gain is the change between the two.
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("E16")
B = base_eq(); ob = outcomes(B); zeta = B.params.zeta_Rl

phis = [0.25, 0.3, 0.35, 0.425, 0.5, 0.6, 0.7]
res = cached("phi_a_sweep") do
    map(phis) do phi
        own = phi == B.params.phi_a
        b = own ? ob : outcomes(ge((phi_a = phi, phi_af = phi); bracket = (0.8, 1.3)))
        r = own ? outcomes(reform_eq(3)) : outcomes(ge((phi_a = phi, phi_af = phi, zeta_Rl = zeta / 3); bracket = (0.8, 1.3)))
        (b = b, r = r)
    end
end
cols = ["base Y/L vs .425", "CF reliance", "CF share", "liq share", "rate", "reform Y/L %", "reform C %", "reform C main5 %", "reform dCFrel", "reform dliq"]
save_table(dir, "pledgeability.csv", cols,
           ["phi_a $phi" => [pct(x.b["wage"], ob["wage"]), x.b["CF reliance"], x.b["CF share"], 100 * x.b["liq share"], x.b["rate"],
                             pct(x.r["wage"], x.b["wage"]), pct(x.r["C"], x.b["C"]), pct(x.r["C main5"], x.b["C main5"]),
                             x.r["CF reliance"] - x.b["CF reliance"], 100 * (x.r["liq share"] - x.b["liq share"])] for (phi, x) in zip(phis, res)])

p1 = small_plot(title = "Gain from zeta/3, %", xlabel = "phi_a (resale value of capital)")
plot!(p1, phis, [pct(x.r["wage"], x.b["wage"]) for x in res], marker = :circle, label = "Y/L")
plot!(p1, phis, [pct(x.r["C"], x.b["C"]) for x in res], marker = :square, label = "consumption")
vline!(p1, [B.params.phi_a], color = :gray, linestyle = :dot, label = "calibrated")
p2 = small_plot(title = "Before the reform", xlabel = "phi_a (resale value of capital)")
plot!(p2, phis, [x.b["CF reliance"] for x in res], marker = :circle, label = "CF reliance")
plot!(p2, phis, [x.b["liq share"] for x in res], marker = :square, label = "liquidated share of defaults")
savefig(plot(p1, p2, layout = (1, 2), size = (900, 340), dpi = 110, bottom_margin = 5Measures.mm), joinpath(dir, "pledgeability.png"))

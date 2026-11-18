############ E06 and E17: cost reform versus bargaining reform ############
# A reorganization reform can cut the cost of the procedure (zeta) and can shift bargaining power from the lender to
# the owner (theta, the lender's share of the reorganized firm). E06: zeta/1, /2, /3 with theta cut by 0 to 60%, each in
# GE, and the theta cut at which the zeta/3 gain is gone. E17: the other direction, "European" procedures with a
# higher zeta (x2, x3) and a weaker lender position.
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("E06")
B = base_eq(); ob = outcomes(B); zeta, theta = B.params.zeta_Rl, B.params.phi_c

tms = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4]
zds = [1.0, 2.0, 3.0, 0.5, 1/3]                  # zeta divided by these: the last two are zeta x 2 and x 3 (E17)
grid = cached("zeta_theta_grid") do
    Dict((zd, tm) => outcomes(zd == 1.0 && tm == 1.0 ? B : ge((zeta_Rl = zeta / zd, phi_c = theta * tm); bracket = (0.8, 1.3)))
         for zd in zds for tm in tms if zd >= 1.0 || tm in (1.0, 0.7, 0.5))
end

name(zd) = zd == 1.0 ? "zeta" : zd > 1 ? "zeta/$(round(Int, zd))" : "zeta x $(round(Int, 1 / zd))"
save_table(dir, "cost_vs_bargaining.csv", SWEEP_COLS,
           ["$(name(zd)), theta x $tm" => sweep_row(grid[(zd, tm)], ob) for zd in zds[1:3] for tm in tms])
save_table(dir, "european_costs.csv", SWEEP_COLS,
           ["$(name(zd)), theta x $tm" => sweep_row(grid[(zd, tm)], ob) for zd in zds[4:5] for tm in (1.0, 0.7, 0.5)])

# break-even: the theta multiple at which the gain of the zeta cut is zero (linear interpolation on the grid)
function break_even(zd, key)
    g = [pct(grid[(zd, tm)][key], ob[key]) for tm in tms]
    i = findfirst(<(0), g)
    (i === nothing || i == 1) && return NaN
    return tms[i-1] + (tms[i] - tms[i-1]) * g[i-1] / (g[i-1] - g[i])
end
for zd in (2.0, 3.0)
    println(name(zd), ": the gain is zero at theta x ", round(break_even(zd, "wage"), digits = 3), " (Y/L), x ",
            round(break_even(zd, "C"), digits = 3), " (consumption)")
end

p1 = small_plot(title = "Y/L, change from the baseline, %", xlabel = "theta, multiple of the baseline", xflip = true)
for zd in zds[1:3]
    plot!(p1, tms, [pct(grid[(zd, tm)]["wage"], ob["wage"]) for tm in tms], marker = :circle, label = name(zd))
end
hline!(p1, [0], color = :black, label = "")
p2 = small_plot(title = "CF reliance", xlabel = "theta, multiple of the baseline", xflip = true)
for zd in zds[1:3]
    plot!(p2, tms, [grid[(zd, tm)]["CF reliance"] for tm in tms], marker = :circle, label = name(zd))
end
savefig(plot(p1, p2, layout = (1, 2), size = (900, 340), dpi = 110, bottom_margin = 5Measures.mm), joinpath(dir, "cost_vs_bargaining.png"))

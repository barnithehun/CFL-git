############ E07: through which channel does the zeta/3 reform work? ############
# Three steps from the baseline to the reform:
#   1 court only:     zeta/3 in the court's decision and the law of motion, but lenders price loans as in the baseline
#                     (Q_FIXED) and the wage stays at 1. Firms' choices cannot change: only who is liquidated does.
#   2 + loan prices:  lenders reprice, the wage still at 1. Not an equilibrium: the entry value exceeds the entry cost.
#   3 + wage (GE):    the wage rises until free entry holds again. This is the reform of the paper.
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("E07")

reset_params!()
B = solve_eq(); qB = LAST_Q[]
z3 = (zeta_Rl = B.params.zeta_Rl / 3,)
Q_FIXED[] = qB
S1 = solve_eq(z3; wage = 1.0)
Q_FIXED[] = nothing
S2 = solve_eq(z3; wage = 1.0)
S3 = reform_eq(3)
println("court only: policies equal the baseline's: ", S1.SumPol[:, [3, 4, 6, 7]] == B.SumPol[:, [3, 4, 6, 7]],
        "; values equal: ", maximum(abs.(S1.SumPol[:, 18] .- B.SumPol[:, 18])))

ob = outcomes(B)
function row(E)
    o = outcomes(E)
    [pct(E.wage, 1.0), pct(E.c_e, B.c_e), pct(o["C"], ob["C"]), pct(o["C main5"], ob["C main5"]), pct(o["firms"], ob["firms"]),
     pct(o["K"], ob["K"]), o["default rate"], o["liq per 100"], 100 * o["liq share"], o["CF reliance"], o["rate small"]]
end
cols = ["wage %", "entry value %", "C %", "C main5 %", "firms %", "K %", "default rate", "liq per 100", "liq share", "CF reliance", "rate small"]
tab = save_table(dir, "channels.csv", cols, ["baseline" => row(B), "1 court only" => row(S1), "2 + loan prices" => row(S2), "3 + wage (GE)" => row(S3)])

plt = small_plot(title = "Consumption, change from the baseline, %", ylabel = "%")
bar!(plt, 1:3, tab[2:4, "C %"], bar_width = 0.6, label = "", xticks = (1:3, ["court only", "+ loan prices", "+ wage (GE)"]))
hline!(plt, [0], color = :black, label = "")
savefig(plt, joinpath(dir, "channels.png"))

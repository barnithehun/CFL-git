############ E01: welfare and the resource table ############
# Steady-state consumption C = Y - I - Psi (main5) in the baseline and under zeta/2 and zeta/3, and where the change
# in C comes from. Entry costs under both conventions: every potential entrant pays c_e (the code's free entry), or
# only the entrants that stay (main5's Psi).
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("E01")

B, R2, R3 = base_eq(), reform_eq(2), reform_eq(3)
wb, w2, w3 = welfare(B), welfare(R2), welfare(R3)
println("law of motion gaps: ", wb["law of motion gap"], " ", w2["law of motion gap"], " ", w3["law of motion gap"])

keys_res = ["Y", "I: net investment", "Psi: fixed operating costs", "Psi: entry costs", "Psi: capital lost in liquidation",
            "Psi: profits lost in liquidation", "Psi: reorganization, fixed (zeta)", "Psi: reorganization, variable", "Psi", "C",
            "Psi: entry costs, main5", "C, main5 entry", "C, reorganization costs as transfers"]
C0 = wb["C"]
save_table(dir, "resources.csv", ["baseline", "zeta/2", "zeta/3", "z/2 pp of C", "z/3 pp of C"],
           [k => [wb[k], w2[k], w3[k], 100 * (w2[k] - wb[k]) / C0, 100 * (w3[k] - wb[k]) / C0] for k in keys_res])

gain(w, k) = pct(w[k], wb[k])
save_table(dir, "gains.csv", ["zeta/2", "zeta/3"],
           ["Y/L = wage, %" => [pct(R2.wage, B.wage), pct(R3.wage, B.wage)],
            "C, %" => [gain(w2, "C"), gain(w3, "C")],
            "C, main5 entry costs, %" => [gain(w2, "C, main5 entry"), gain(w3, "C, main5 entry")],
            "C, reorganization costs as transfers, %" => [gain(w2, "C, reorganization costs as transfers"), gain(w3, "C, reorganization costs as transfers")],
            "producing firms, %" => [gain(w2, "firms"), gain(w3, "firms")],
            "potential entrants, %" => [gain(w2, "potential entrants"), gain(w3, "potential entrants")],
            "entrants that stay, %" => [gain(w2, "entrants that stay"), gain(w3, "entrants that stay")],
            "K, %" => [gain(w2, "K"), gain(w3, "K")]])

# figure: contributions to the change in C, in percent of baseline C (costs enter with a minus sign)
items = ["Output" => "Y", "Investment" => "I: net investment", "Fixed op. costs" => "Psi: fixed operating costs",
         "Entry costs" => "Psi: entry costs", "Liquidation losses" => ["Psi: capital lost in liquidation", "Psi: profits lost in liquidation"],
         "Reorg. cost zeta" => "Psi: reorganization, fixed (zeta)", "Reorg. variable" => "Psi: reorganization, variable"]
contrib(w, k) = (k isa String ? w[k] - wb[k] : sum(w[j] - wb[j] for j in k)) * (k == "Y" ? 1 : -1) * 100 / C0
labels = [first(i) for i in items]
c2, c3 = [contrib(w2, last(i)) for i in items], [contrib(w3, last(i)) for i in items]
plt = small_plot(title = "Change in consumption, % of baseline C", ylabel = "pp", xrotation = 30, bottom_margin = 8Measures.mm)
bar!(plt, (1:length(items)) .- 0.2, c2, bar_width = 0.38, label = "zeta/2 (total $(round(sum(c2), digits = 2)))", xticks = (1:length(items), labels))
bar!(plt, (1:length(items)) .+ 0.2, c3, bar_width = 0.38, label = "zeta/3 (total $(round(sum(c3), digits = 2)))")
hline!(plt, [0], color = :black, label = "")
savefig(plt, joinpath(dir, "consumption_change.png"))

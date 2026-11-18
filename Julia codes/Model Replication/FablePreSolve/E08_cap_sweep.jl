############ E08 and E13: who is eligible for the cheaper procedure? ############
# E08: zeta/3 only for borrowers with debt below a cap, for caps from the smallest borrowers to everyone, each in GE.
# The caps sit between two points of the debt grid, so no firm is exactly at the cap. Bunching: the mass of borrowers
# at the grid point just below the cap, against the baseline's mass there.
# E13: the same cut for the same share of borrowers, chosen by debt or by capital (Kim 2021: the Chapter 11 cost
# disadvantage is in small-asset firms, not small-debt firms).
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("E08")
B = base_eq(); ob = outcomes(B); zeta = B.params.zeta_Rl
S, c = B.SumPol, cont_mass(B)
b, k = S[1:end-1, 4], S[1:end-1, 3]
bor = b .> 0

# caps at the geometric midpoints of the debt grid, where the cumulative share of baseline borrowers is nearest to the targets
vals = sort(unique(b[bor .& (c .> 0)]))
cum = cumsum([sum(c[b .== v]) for v in vals]) ./ sum(c[bor])
targets = [0.1, 0.2, 0.35, 0.5, 0.68, 0.8, 0.9, 0.97]
idx = unique([clamp(argmin(abs.(cum .- t)), 1, length(vals) - 1) for t in targets])
caps = [sqrt(vals[i] * vals[i+1]) for i in idx]

borrower_share(E, v) = (cE = cont_mass(E); bE = E.SumPol[1:end-1, 4]; sum(cE[bE .== v]) / sum(cE[bE .> 0]))
res = cached("cap_sweep") do
    map(zip(idx, caps)) do (i, cap)
        E = ge((zeta_Rs = zeta / 3, Fcut = cap))
        (o = outcomes(E), below = borrower_share(E, vals[i]), above = borrower_share(E, vals[i+1]))
    end
end
rows = ["cap $(round(Int, cap)) ($(round(Int, 100 * cum[i]))% eligible)" =>
        [sweep_row(r.o, ob)[1:4]; 100 * cum[i]; 100 * r.below; 100 * borrower_share(B, vals[i]); 100 * r.above; 100 * borrower_share(B, vals[i+1])]
        for (i, cap, r) in zip(idx, caps, res)]
oU = outcomes(reform_eq(3))
push!(rows, "no cap (all firms)" => [sweep_row(oU, ob)[1:4]; 100.0; NaN; NaN; NaN; NaN])
save_table(dir, "cap_sweep.csv", ["Y/L %", "C %", "C main5 %", "firms %", "eligible %", "below: reform", "below: base", "above: reform", "above: base"], rows)

p1 = small_plot(title = "Gain by eligibility cap, %", xlabel = "eligible share of baseline borrowers, %")
plot!(p1, [100 * cum[idx]; 100], [[pct(r.o["wage"], ob["wage"]) for r in res]; pct(oU["wage"], ob["wage"])], marker = :circle, label = "Y/L")
plot!(p1, [100 * cum[idx]; 100], [[pct(r.o["C"], ob["C"]) for r in res]; pct(oU["C"], ob["C"])], marker = :square, label = "consumption")
p2 = small_plot(title = "Borrowers just below the cap, % of borrowers", xlabel = "eligible share of baseline borrowers, %")
plot!(p2, 100 * cum[idx], [100 * borrower_share(B, vals[i]) for i in idx], marker = :circle, label = "baseline")
plot!(p2, 100 * cum[idx], [100 * r.below for r in res], marker = :square, label = "with the cap")
savefig(plot(p1, p2, layout = (1, 2), size = (900, 340), dpi = 110, bottom_margin = 5Measures.mm), joinpath(dir, "cap_sweep.png"))

# E13: eligibility by debt (the model's Fcut) or by capital, for the same share of baseline borrowers
share_debt = sum(c[bor .& (b .< B.params.Fcut)]) / sum(c[bor])
kvals = sort(unique(k[bor .& (c .> 0)]))
kcum = cumsum([sum(c[bor .& (k .== v)]) for v in kvals]) ./ sum(c[bor])
j = clamp(argmin(abs.(kcum .- share_debt)), 1, length(kvals) - 1)
Fcut_k = sqrt(kvals[j] * kvals[j+1])
println("\nE13: eligible share of borrowers by debt ", round(share_debt, digits = 3), " (Fcut ", B.params.Fcut, "); by capital ",
        round(kcum[j], digits = 3), " (k' < ", round(Fcut_k, digits = 1), ")")
elig = cached("eligibility") do
    (debt = outcomes(ge((zeta_Rs = zeta / 3,))), capital = outcomes(ge((zeta_Rs = zeta / 3, size_by_debt = 0.0, Fcut = Fcut_k))))
end
extra(o) = [o["reorganizations"] - ob["reorganizations"], o["reo per 100"], o["liq per 100"]]
save_table(dir, "eligibility.csv", [SWEEP_COLS; "added reorg."; "reo per 100"; "liq per 100"],
           ["baseline" => [sweep_row(ob, ob); extra(ob)], "by debt" => [sweep_row(elig.debt, ob); extra(elig.debt)],
            "by capital" => [sweep_row(elig.capital, ob); extra(elig.capital)], "all firms" => [sweep_row(oU, ob); extra(oU)]])

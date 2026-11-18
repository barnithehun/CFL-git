############ E09: how often is the court wrong? ############
# The court sees the reorganization value with noise (sigma_v). Without noise it would liquidate exactly when
# L + zeta >= phi_c_hh * V (L = phi_a (1 - delta) k'). A mistaken liquidation is a liquidation of a firm the noiseless
# court would reorganize; a mistaken reorganization is the reverse. The value lost in a mistake is the gap between
# the two sides. Compare: Antill (2022) finds many liquidations that reduce recoveries; Dou et al. (2021) find excess
# liquidation and continuation to be small; Jenkins and Smith (2014) put the loss at up to 0.28% of going-concern value.
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("E09")

function mistakes(E, group)
    ev = events(E); c = cont_mass(E); p = E.params
    w = c[ev.s] .* ev.prob
    keep = (ev.b .> 0) .& (group == :all ? trues(length(w)) : group == :small ? .!ev.large : ev.large)
    d = (w .* ev.pdef)[keep]
    liq, k, v, z = ev.liq[keep], ev.k[keep], ev.vreo[keep], ev.zeta[keep]
    L = p.phi_a .* (1 - p.delta) .* k
    gc = p.phi_c_hh .* v .- z                          # what reorganization is worth to lender and owner
    should = L .>= gc                                  # the noiseless court liquidates
    mliq = d .* liq .* .!should
    mreo = d .* (1 .- liq) .* should
    value = sum(d .* max.(L, gc))                      # value of the defaulting firms under the right decisions
    C = welfare(E)["C"]; sc = 10000 / (statevecs(E, ev).L' * c)
    return ["defaults per 100 borrowers of the group" => 100 * sum(d) / sum(c[(E.SumPol[1:ev.n, 4] .> 0) .& (group == :all ? trues(ev.n) :
                                                          (((p.size_by_debt == 1 ? E.SumPol[1:ev.n, 4] : E.SumPol[1:ev.n, 3]) .>= p.Fcut) .== (group == :large)))]),
            "liquidated share of defaults, %" => 100 * sum(d .* liq) / sum(d),
            "  the noiseless court's share, %" => 100 * sum(d[should]) / sum(d),
            "mistaken liquidations, % of liquidations" => 100 * sum(mliq) / sum(d .* liq),
            "mistaken reorganizations, % of reorganizations" => 100 * sum(mreo) / sum(d .* (1 .- liq)),
            "all mistakes, % of defaults" => 100 * (sum(mliq) + sum(mreo)) / sum(d),
            "value lost in mistaken liquidations, % of value" => 100 * sum(mliq .* (gc .- L)) / value,
            "value lost in mistaken reorganizations, % of value" => 100 * sum(mreo .* (L .- gc)) / value,
            "value lost in all mistakes, % of consumption" => 100 * sc * (sum(mliq .* (gc .- L)) + sum(mreo .* (L .- gc))) / C]
end

B, R3 = base_eq(), reform_eq(3)
cols = [mistakes(B, :all), mistakes(B, :small), mistakes(B, :large), mistakes(R3, :all), mistakes(R3, :small), mistakes(R3, :large)]
tab = save_table(dir, "mistakes.csv", ["base all", "base small", "base large", "z/3 all", "z/3 small", "z/3 large"],
                 [first(cols[1][i]) => [last(col[i]) for col in cols] for i in eachindex(cols[1])])

# figure: defaults split into right and wrong decisions, small and large borrowers, baseline
function split4(col)
    liqsh, ml, mr = last(col[2]) / 100, last(col[4]) / 100, last(col[5]) / 100
    100 .* [liqsh * (1 - ml), liqsh * ml, (1 - liqsh) * mr, (1 - liqsh) * (1 - mr)]
end
parts = hcat(split4(cols[2]), split4(cols[3]), split4(cols[5]), split4(cols[6]))'
plt = small_plot(title = "How defaults are resolved, % of defaults", ylims = (0, 100), legend = :outerright, size = (720, 360))
acc = zeros(4)
for (j, lab) in enumerate(["liquidated, right", "liquidated, mistake", "reorganized, mistake", "reorganized, right"])
    bar!(plt, 1:4, acc .+ parts[:, j], fillrange = acc, bar_width = 0.6, label = lab,
         xticks = (1:4, ["small, base", "large, base", "small, z/3", "large, z/3"]))
    acc .+= parts[:, j]
end
savefig(plt, joinpath(dir, "mistakes.png"))

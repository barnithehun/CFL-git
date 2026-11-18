############ E04 and E15: the size gradient ############
# E04: borrowers in bins of debt b' (the bins are cut at the deciles of the baseline's borrowers, so both economies are
# compared at the same debt levels): liquidated share of defaults, CF reliance, interest rate, default rate.
# E15: all firms in deciles of employment: share with debt and CF reliance (is CF reliance U-shaped in size?).
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("E04")
B, R3 = base_eq(), reform_eq(3)

# bin edges: the b' values at which the cumulative mass of baseline borrowers passes each decile
function debt_edges(E; q = 0.1:0.1:0.9)
    c, b = cont_mass(E), E.SumPol[1:end-1, 4]
    vals = sort(unique(b[(b .> 0) .& (c .> 0)]))
    cum = cumsum([sum(c[b .== v]) for v in vals]) ./ sum(c[b .> 0])
    return unique([vals[findfirst(>=(x), cum)] for x in q])
end
edges = debt_edges(B)
bin_of(b) = searchsortedfirst(edges, b)             # 1 .. length(edges) + 1; a bin holds b' up to and including its edge

function by_debt(E)
    ev = events(E); c = cont_mass(E); S = E.SumPol
    w = c[ev.s] .* ev.prob
    out = Vector{Vector{Float64}}()
    for j in 1:length(edges) + 1
        st = (S[1:ev.n, 4] .> 0) .& (bin_of.(S[1:ev.n, 4]) .== j)
        rows = st[ev.s]
        d = w[rows] .* ev.pdef[rows]
        ss = state_stats(E, c .* st)
        push!(out, [100 * sum(c[st]) / sum(c[S[1:ev.n, 4] .> 0]), 100 * sum(d .* ev.liq[rows]) / sum(d), ss.tau, ss.rate,
                    100 * sum(d) / sum(c[st]), ss.d2c])
    end
    return out
end
b0, b3 = by_debt(B), by_debt(R3)
labels = [j <= length(edges) ? "b' <= $(round(Int, edges[j]))" : "b' > $(round(Int, edges[end]))" for j in 1:length(edges) + 1]
save_table(dir, "by_debt.csv", ["share base", "share z/3", "liq base", "liq z/3", "CFrel base", "CFrel z/3", "rate base", "rate z/3", "def base", "def z/3"],
           [labels[j] => [b0[j][1], b3[j][1], b0[j][2], b3[j][2], b0[j][3], b3[j][3], b0[j][4], b3[j][4], b0[j][5], b3[j][5]] for j in eachindex(labels)])
println("Fcut (small/large split) = ", B.params.Fcut)

x = 1:length(labels)
p1 = small_plot(title = "Liquidated share of defaults, %", xlabel = "debt bin (baseline deciles)")
plot!(p1, x, [v[2] for v in b0], marker = :circle, label = "baseline"); plot!(p1, x, [v[2] for v in b3], marker = :square, label = "zeta/3")
p2 = small_plot(title = "CF reliance", xlabel = "debt bin (baseline deciles)")
plot!(p2, x, [v[3] for v in b0], marker = :circle, label = "baseline"); plot!(p2, x, [v[3] for v in b3], marker = :square, label = "zeta/3")
savefig(plot(p1, p2, layout = (1, 2), size = (900, 340), dpi = 110, bottom_margin = 5Measures.mm), joinpath(dir, "by_debt.png"))

# E15: deciles of employment over all continuing firms
function by_employment(E)
    c, S = cont_mass(E), E.SumPol
    lab = S[1:end-1, 10]
    ord = sortperm(lab)
    cum = cumsum(c[ord]) ./ sum(c)
    dec = zeros(Int, length(c)); dec[ord] = min.(10, floor.(Int, cum .* 10 .- 1e-9) .+ 1)
    return [(st = dec .== j; ss = state_stats(E, c .* st); [ss.lab, 100 * ss.borrowing, ss.tau,
             sum(c[st] .* S[1:end-1, 4][st] .* S[1:end-1, 17][st]) / max(sum(c[st] .* S[1:end-1, 4][st]), 1e-12), ss.d2c]) for j in 1:10]
end
e0, e3 = by_employment(B), by_employment(R3)
save_table(dir, "by_employment.csv", ["employees", "borrow % base", "borrow % z/3", "CFrel base", "CFrel z/3", "CFshare base", "CFshare z/3"],
           ["decile $j" => [e0[j][1], e0[j][2], e3[j][2], e0[j][3], e3[j][3], e0[j][4], e3[j][4]] for j in 1:10])
p3 = small_plot(title = "CF reliance of borrowers", xlabel = "employment decile")
plot!(p3, 1:10, [v[3] for v in e0], marker = :circle, label = "baseline"); plot!(p3, 1:10, [v[3] for v in e3], marker = :square, label = "zeta/3")
p4 = small_plot(title = "Firms with debt, %", xlabel = "employment decile")
plot!(p4, 1:10, [v[2] for v in e0], marker = :circle, label = "baseline"); plot!(p4, 1:10, [v[2] for v in e3], marker = :square, label = "zeta/3")
savefig(plot(p3, p4, layout = (1, 2), size = (900, 340), dpi = 110, bottom_margin = 5Measures.mm), joinpath(dir, "by_employment.png"))

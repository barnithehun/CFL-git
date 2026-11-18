############ E14: which firms does the reform save? ############
# In the zeta/3 economy, every defaulting firm has a liquidation probability under the new cost and one under the old
# cost (same firm, same value). The difference is the mass of "saved" firms: reorganized now, liquidated before.
# Are they weak firms kept alive (zombies)? Compare their productivity and capital, and their next five years, with
# the firms reorganized under both costs, the firms still liquidated, and the average firm.
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("E14")
B, R = base_eq(), reform_eq(3)
ev = events(R); sv = statevecs(R, ev); K = kernel(ev); c = cont_mass(R)
w = c[ev.s] .* ev.prob
bor = ev.b .> 0
liq_old = liq_prob(R, ev.vreo, ev.k, B.params.zeta_Rl)
d = w .* ev.pdef .* bor
groups = ["saved by the reform" => d .* (liq_old .- ev.liq), "reorganized under both" => d .* (1 .- liq_old), "still liquidated" => d .* ev.liq]
avg = state_stats(R, c)
e_next = ev.e_vals[ev.en]
ent = entrant_dist(R)

# where a group of reorganized firms lands (x' + b')
function landing_dist(g)
    d0 = zeros(ev.n)
    for i in eachindex(g)
        g[i] == 0 && continue
        d0[ev.r1[i]] += g[i] * ev.rw[i]
        ev.r2[i] != 0 && (d0[ev.r2[i]] += g[i] * (1 - ev.rw[i]))
    end
    return d0
end
rows = Pair{String, Vector{Float64}}[]
paths = Dict{String, Any}()
for (name, g) in groups
    f = name == "still liquidated" ? nothing : follow(R, K, sv, landing_dist(g); T = 5)
    paths[name] = f
    wm(x) = wmean(x, g, g .> 0)
    push!(rows, name => [100 * sum(g) / sum(d), wm(e_next) / avg.e, wm(ev.k) / avg.k, wm(ev.b ./ max.(ev.k, 1e-12)), wm(ev.tau),
                         100 * wm(Float64.(ev.large)), f === nothing ? NaN : 100 * f.alive[4], f === nothing ? NaN : 100 * f.cdef[4],
                         f === nothing ? NaN : 100 * f.alive[6]])
end
fa = follow(R, K, sv, K.T * c; T = 5); paths["average firm"] = fa
push!(rows, "average firm" => [NaN, 1.0, 1.0, avg.d2c, avg.tau, NaN, 100 * fa.alive[4], 100 * fa.cdef[4], 100 * fa.alive[6]])
fe = follow(R, K, sv, K.T * (ent ./ sum(ent)); T = 5); paths["entrants"] = fe
push!(rows, "entrants" => [NaN, sum(ent .* R.SumPol[1:ev.n, 2]) / sum(ent) / avg.e, sum(ent .* R.SumPol[1:ev.n, 3]) / sum(ent) / avg.k, NaN, NaN, NaN,
                           100 * fe.alive[4], 100 * fe.cdef[4], 100 * fe.alive[6]])
save_table(dir, "saved_firms.csv", ["% of defaults", "productivity", "capital", "debt/capital", "CF reliance", "large %", "alive 3y %", "default 3y %", "alive 5y %"], rows)

plt = small_plot(title = "Still operating, % (zeta/3 economy)", xlabel = "years after the default (or from now)")
for (name, sty) in ["saved by the reform" => :solid, "reorganized under both" => :dash, "average firm" => :dot, "entrants" => :dashdot]
    plot!(plt, 0:5, 100 .* paths[name].alive, linestyle = sty, marker = :circle, markersize = 3, label = name)
end
savefig(plt, joinpath(dir, "saved_firms.png"))

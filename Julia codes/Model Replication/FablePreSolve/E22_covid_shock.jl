############ E22: a one-off revenue shock with and without the reform ############
# In year 0 every firm unexpectedly loses the share s of its revenue (wages are paid in full), as in a lockdown. The
# economy is in its steady state before: the baseline, or the economy that already has the zeta/3 law. Firms arrive
# with less cash; those that end up at a point where they default are liquidated or reorganized under the law in
# force. The shock year is computed exactly from the model's events.
# The recovery after year 0 uses E21's approximation (path_static: steady-state policies, steady-state entry, the wage
# clears the labor market given the capital in place), because the constant-wage path does not exist (see E21).
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("E22")
T = 60

function shock_path(E, s)
    ev = events(E); sv = statevecs(E, ev); K = kernel(ev)
    evS = events(E; rev_shock = s); svS = statevecs(E, evS); KS = kernel(evS)
    c0 = cont_mass(E); sc = 10000 / (sv.L' * c0); c = sc .* c0
    w = welfare(E; ev = ev, sv = sv)
    m_ss, nu, beta = w["potential entrants"], E.params.nu, E.params.beta
    P = path_static(E, KS.T * c, m_ss; T = T + 1, ev = ev, sv = sv, K = K)
    r0 = resources(E, svS, c, P.c[1], m_ss)
    C = [r0["C"] - s * r0["Y"]]                          # the shock year: revenue s * Y is lost
    wage = [E.wage]
    for t in 1:T
        r = resources(E, sv, P.c[t], P.c[t+1], m_ss)
        push!(C, r["C"] - r["Y"] + P.wage[t] * 10000 / nu); push!(wage, P.wage[t])
    end
    return (C = C, Css = w["C"], wage = wage, firms = [sum(c); [sum(P.c[t]) for t in 1:T]], firms_ss = w["firms"],
            liq0 = svS.liq' * c, liq_ss = sv.liq' * c, reo0 = svS.reo' * c, reo_ss = sv.reo' * c, def0 = svS.def' * c, def_ss = sv.def' * c,
            liqK0 = svS.liqK' * c, liqK_ss = sv.liqK' * c, K = sv.K' * c, K_after = sv.K' * P.c[1],
            pv = sum(beta .^ (0:T) .* (C .- w["C"])) / w["C"])  # discounted consumption loss, in years of consumption
end

B, R = base_eq(), reform_eq(3)
shocks = [0.1, 0.2, 0.3, 0.5]
res = Dict((lab, s) => shock_path(E, s) for (lab, E) in ["baseline" => B, "zeta/3" => R] for s in shocks)
rows = Pair{String, Vector{Float64}}[]
for s in shocks, lab in ["baseline", "zeta/3"]
    x = res[(lab, s)]
    push!(rows, "shock $(round(Int, 100s))%, $lab" =>
          [100 * x.def0 / x.firms_ss, 100 * x.def_ss / x.firms_ss, 100 * x.liq0 / x.firms_ss, 100 * x.liq_ss / x.firms_ss, 100 * x.reo0 / x.firms_ss,
           100 * x.liqK0 / x.K, 100 * x.liqK_ss / x.K, pct(x.K_after, x.K), pct(x.wage[2], x.wage[1]), pct(minimum(x.firms), x.firms_ss)])
end
save_table(dir, "shock.csv", ["defaults %", "(normal)", "liquidated %", "(normal)", "reorganized %", "capital liq. %", "(normal)", "capital chosen %",
                              "wage year 1 %", "firms trough %"], rows)

s = 0.3
p1 = small_plot(title = "Liquidated in the shock year, per 100 firms", xlabel = "revenue lost, %")
p2 = small_plot(title = "Capital liquidated in the shock year, % of capital", xlabel = "revenue lost, %")
for (lab, sty) in ["baseline" => :solid, "zeta/3" => :dash]
    plot!(p1, 100 .* shocks, [100 * (res[(lab, x)].liq0 - res[(lab, x)].liq_ss) / res[(lab, x)].firms_ss for x in shocks], linestyle = sty, marker = :circle,
          label = "$lab: extra liquidations")
    plot!(p1, 100 .* shocks, [100 * (res[(lab, x)].reo0 - res[(lab, x)].reo_ss) / res[(lab, x)].firms_ss for x in shocks], linestyle = sty, marker = :square,
          label = "$lab: extra reorganizations")
    plot!(p2, 100 .* shocks, [100 * res[(lab, x)].liqK0 / res[(lab, x)].K for x in shocks], linestyle = sty, marker = :circle, label = lab)
end
savefig(plot(p1, p2, layout = (1, 2), size = (900, 340), dpi = 110, bottom_margin = 5Measures.mm), joinpath(dir, "shock.png"))

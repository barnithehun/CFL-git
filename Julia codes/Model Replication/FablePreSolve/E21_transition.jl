############ E21: the transition after the reform ############
# The reform (zeta/2 or zeta/3) is announced without warning at the end of period 0's bankruptcies: the firms arrive in
# period 0 as in the baseline steady state, and from their choices in period 0 on the new law holds.
#
# Step 1, the constant-wage path. With free entry and a fixed interest rate, the wage is pinned down by the entry
# condition alone as long as some firms enter. If that held in every period, the wage would jump to its new level at
# once and only the firm distribution would move (path()). It does not hold: entrants are small and grow fast, so the
# entrants needed to clear the labor market in one year demand too much labor in the next years, and entry hits zero.
# The true path therefore has a moving wage, which needs value functions that change over time: a new solver.
#
# Step 2, an approximation (path_static): the firms follow the new steady state's policies, entry is at its new
# steady-state level from period 0 on, and each period the wage clears the labor market given the capital in place.
# It shows how fast the firm distribution turns over; firms do not foresee the wage deviations, so it is not exact.
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("E21")
B = base_eq()
evB = events(B); svB = statevecs(B, evB); KB = kernel(evB)
cB = cont_mass(B); scB = 10000 / (svB.L' * cB)
CB = welfare(B; ev = evB, sv = svB)["C"]
beta, nu, T = B.params.beta, B.params.nu, 120

results = Dict{Int, Any}()
for dv in (2, 3)
    R = reform_eq(dv)
    evR = events(R); svR = statevecs(R, evR); KR = kernel(evR)
    wR = welfare(R; ev = evR, sv = svR)
    d0 = scB .* (KB.T * cB)                           # incumbents arriving in period 0
    m_ss = wR["potential entrants"]

    # step 1: the constant-wage path and its entry
    P1 = path(R, d0; T = 30, ev = evR, sv = svR, K = KR)
    println("\nzeta/$dv, constant-wage path: years with zero entry among the first 30: ", P1.corners,
            "; entry relative to the new steady state, years 0-10: ", round.(P1.m[1:11] ./ m_ss, digits = 2))

    # step 2: the approximation
    P = path_static(R, d0, m_ss; T = T + 1, ev = evR, sv = svR, K = KR)
    # period 0: the baseline's firms produce at the baseline wage; period t >= 1: the firms P.c[t] produce at wage P.wage[t]
    C = [resources(B, svB, scB .* cB, P.c[1], m_ss; sv_next = svR)["C"]]
    wage = [B.wage]
    for t in 1:T
        r = resources(R, svR, P.c[t], P.c[t+1], m_ss)
        push!(C, r["C"] - r["Y"] + P.wage[t] * 10000 / nu); push!(wage, P.wage[t])
    end
    disc = beta .^ (0:T)
    pv = (sum(disc .* (C .- CB)) + beta^(T + 1) / (1 - beta) * (wR["C"] - CB)) / (CB / (1 - beta))
    firms = [scB * sum(cB); [sum(P.c[t]) for t in 1:T]]
    closed(x, x0, x1) = 100 .* (x .- x0) ./ (x1 - x0)
    half(x, x0, x1) = (i = findfirst(>=(50), closed(x, x0, x1)); i === nothing ? NaN : i - 1)
    results[dv] = (C = C, wage = wage, firms = firms, pv = pv, ss = pct(wR["C"], CB), wage_ss = R.wage, firms_ss = wR["firms"],
                   half_wage = half(wage, B.wage, R.wage), half_firms = half(firms, firms[1], wR["firms"]), corners = P1.corners)
    println("zeta/$dv, approximation: wage gap closed after 1, 2, 5, 10, 20 years, %: ", round.(closed(wage, B.wage, R.wage)[[2, 3, 6, 11, 21]], digits = 1),
            "; firm mass gap closed: ", round.(closed(firms, firms[1], wR["firms"])[[2, 3, 6, 11, 21]], digits = 1))
    println("   consumption, % from the baseline, years 0-10: ", round.(pct.(C[1:11], CB), digits = 2), "; year $T against the new steady state: ",
            round(pct(C[end], wR["C"]), digits = 4), "%")
    println("   steady-state consumption gain ", round(pct(wR["C"], CB), digits = 3), "%, with the transition (discounted) ", round(100 * pv, digits = 3), "%")
end

df = DataFrame(year = 0:T)
for dv in (2, 3)
    x = results[dv]
    df[!, "C % from baseline, zeta/$dv"] = pct.(x.C, CB); df[!, "wage % from baseline, zeta/$dv"] = pct.(x.wage, B.wage)
    df[!, "firms, zeta/$dv"] = x.firms
end
CSV.write(joinpath(dir, "transition.csv"), df)
save_table(dir, "transition_welfare.csv", ["zeta/2", "zeta/3"],
           ["constant-wage path: years with zero entry (of 30)" => [results[2].corners, results[3].corners],
            "steady-state wage (Y/L) gain, %" => [pct(results[2].wage_ss, B.wage), pct(results[3].wage_ss, B.wage)],
            "years until half of the wage gain" => [results[2].half_wage, results[3].half_wage],
            "years until half of the rise in firms" => [results[2].half_firms, results[3].half_firms],
            "steady-state consumption gain, %" => [results[2].ss, results[3].ss],
            "consumption gain with the transition, %" => [100 * results[2].pv, 100 * results[3].pv],
            "consumption in year 0, % from baseline" => [pct(results[2].C[1], CB), pct(results[3].C[1], CB)],
            "consumption in year 5, % from baseline" => [pct(results[2].C[6], CB), pct(results[3].C[6], CB)]])

yrs = 0:40
p1 = small_plot(title = "Wage = Y/L, % from the baseline (approximation)", xlabel = "years after the reform")
p2 = small_plot(title = "Consumption, % from the baseline (approximation)", xlabel = "years after the reform")
for dv in (2, 3)
    plot!(p1, yrs, pct.(results[dv].wage[1:41], B.wage), label = "zeta/$dv")
    hline!(p1, [pct(results[dv].wage_ss, B.wage)], linestyle = :dot, color = :gray, label = dv == 2 ? "new steady states" : "")
    plot!(p2, yrs, pct.(results[dv].C[1:41], CB), label = "zeta/$dv")
    hline!(p2, [results[dv].ss], linestyle = :dot, color = :gray, label = dv == 2 ? "new steady states" : "")
end
savefig(plot(p1, p2, layout = (1, 2), size = (900, 340), dpi = 110, bottom_margin = 5Measures.mm), joinpath(dir, "transition.png"))

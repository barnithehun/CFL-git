############ E11: firms around a default ############
# Year 0 is the year of the default. Years -5 to -1: the firms that will be reorganized (or liquidated) in year 0,
# at their choices 5 to 1 years before (the last choice before the default is year -1). Years 0 to 5: the reorganized
# firms after their fresh start. Averages are relative to all continuing firms of the same economy.
# Before the event the weights are exact: c(s) * P(the event happens j years later | state s).
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("E11")

function event_study(E; J = 5)
    ev = events(E); sv = statevecs(E, ev); K = kernel(ev); cont = continuing(E); c = cont_mass(E)
    avg = state_stats(E, c)
    stat(wts) = (st = state_stats(E, wts); [st.k / avg.k, st.e / avg.e, st.d2c, st.tau, st.rate, 100 * st.borrowing])
    pre = Dict{Symbol, Vector{Vector{Float64}}}()
    for (ename, h1) in [:reorganized => sv.reo, :liquidated => sv.liq]
        h, out = copy(h1), Vector{Vector{Float64}}()
        for j in 1:J
            push!(out, stat(c .* h))
            h = K.T' * (cont .* h)                    # probability that the event happens one more year ahead
        end
        pre[ename] = reverse(out)                     # years -J .. -1
    end
    post, alive = Vector{Vector{Float64}}(), Float64[]
    d = K.Treo * c; n0 = sum(d)
    for t in 0:J
        cc = d .* cont
        push!(post, stat(cc)); push!(alive, 100 * sum(cc) / n0)
        d = K.T * cc
    end
    return (pre = pre, post = post, alive = alive)
end
names_ = ["capital (rel.)", "productivity (rel.)", "debt/capital", "CF reliance", "rate %", "with debt %"]
e0, e3 = event_study(base_eq()), event_study(reform_eq(3))

df = DataFrame(year = -5:5)
for (lab, e) in ["base" => e0, "z/3" => e3], (j, nm) in enumerate(names_)
    df[!, "$nm, reorganized, $lab"] = [[r[j] for r in e.pre[:reorganized]]; [r[j] for r in e.post]]
    df[!, "$nm, liquidated, $lab"] = [[r[j] for r in e.pre[:liquidated]]; fill(NaN, 6)]
end
df[!, "reorganized still continuing %, base"] = [fill(NaN, 5); e0.alive]
df[!, "reorganized still continuing %, z/3"] = [fill(NaN, 5); e3.alive]
CSV.write(joinpath(dir, "event_study.csv"), df)
println("\nyear | capital, productivity, debt/capital, CF reliance: reorganized base, reorganized z/3, liquidated base")
for (i, y) in enumerate(-5:5)
    println(lpad(y, 4), " | ", join([string(round(df[i, "$nm, reorganized, base"], digits = 2), " ", round(df[i, "$nm, reorganized, z/3"], digits = 2), " ",
                                     round(df[i, "$nm, liquidated, base"], digits = 2)) for nm in names_[1:4]], " | "))
end
println("reorganized firms still continuing, %: base ", round.(e0.alive, digits = 1), " z/3 ", round.(e3.alive, digits = 1))

ps = []
for (j, nm) in enumerate(names_[1:4])
    p = small_plot(title = nm, xlabel = "years from the default")
    plot!(p, -5:5, df[!, "$nm, reorganized, base"], marker = :circle, markersize = 3, label = "reorganized, baseline")
    plot!(p, -5:5, df[!, "$nm, reorganized, z/3"], marker = :square, markersize = 3, linestyle = :dash, label = "reorganized, zeta/3")
    plot!(p, -5:-1, df[1:5, "$nm, liquidated, base"], marker = :diamond, markersize = 3, linestyle = :dot, label = "liquidated, baseline")
    vline!(p, [-0.5], color = :gray, label = "")
    push!(ps, p)
end
savefig(plot(ps..., layout = (2, 2), size = (900, 600), dpi = 110, bottom_margin = 4Measures.mm), joinpath(dir, "event_study.png"))

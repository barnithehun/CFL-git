############ E10: the life cycle of an entry cohort ############
# One cohort of entrants followed through the law of motion for 30 years: survival, size, borrowing, CF reliance,
# interest rate and default by age, in the baseline and under zeta/3 (each at its own GE wage).
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("E10")

function cohort(E; T = 30)
    ev = events(E); sv = statevecs(E, ev); K = kernel(ev); cont = continuing(E)
    c = entrant_dist(E); c = c ./ sum(c)             # entrants that stay, mass 1; they decide in the entry period (age 0)
    rows = Vector{Vector{Float64}}()
    for age in 0:T
        st = state_stats(E, c)
        push!(rows, [st.mass, st.lab, st.k, 100 * st.borrowing, st.tau, st.rate, st.d2c,
                     100 * (sv.def' * c) / st.mass, 100 * (sv.liq' * c) / st.mass, st.e])
        c = (K.T * c) .* cont
    end
    return rows
end
names_ = ["surviving", "employees", "capital", "with debt %", "CF reliance", "rate %", "debt/capital", "default %", "liquidated %", "productivity"]
c0, c3 = cohort(base_eq()), cohort(reform_eq(3))
df = DataFrame(age = 0:30)
for (j, nm) in enumerate(names_)
    df[!, nm * ", base"] = [r[j] for r in c0]; df[!, nm * ", z/3"] = [r[j] for r in c3]
end
CSV.write(joinpath(dir, "age_profiles.csv"), df)
println("\nage profiles (base | zeta/3)")
for a in [0, 1, 2, 3, 5, 10, 20, 30]
    println(rpad("age $a", 8), join([lpad(string(round(c0[a+1][j], digits = 2), " | ", round(c3[a+1][j], digits = 2)), 22) for j in [1, 2, 4, 5, 6, 8]]))
end
println(rpad("", 8), join(lpad.(names_[[1, 2, 4, 5, 6, 8]], 22)))

ps = []
for j in [1, 2, 4, 5, 6, 8]
    p = small_plot(title = names_[j], xlabel = "age")
    plot!(p, 0:30, [r[j] for r in c0], label = "baseline"); plot!(p, 0:30, [r[j] for r in c3], label = "zeta/3", linestyle = :dash)
    push!(ps, p)
end
savefig(plot(ps..., layout = (2, 3), size = (1000, 560), dpi = 110, bottom_margin = 4Measures.mm), joinpath(dir, "age_profiles.png"))

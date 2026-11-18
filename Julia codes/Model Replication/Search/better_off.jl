# better_off.jl <run id> [name=value ...]: are small firms better off after the zeta/2 reform? At a run's calibrated values
# (runs/<id>/CALIB.csv and the external settings in its DIAG.csv), optionally on another grid (e.g. k_size=80 b_size=80),
# solves the baseline and the zeta/2 reform in GE and writes better_off/<id>_<grid>.csv with
#   - the fit at these values on this grid (five aggregates, small and large moments; Fcut refitted to the data's small share)
#   - small firms before and after: averages over small firms in each equilibrium ("total"), and the reform policies at the
#     baseline's small-firm states weighted by the baseline distribution ("same states": incumbents, no entry or exit effect)
#   - the mean firm value V of small firms (SumPol column 18) in both measures
#   - "p": the same-state measure with the reform at the baseline wage (partial equilibrium: the credit effect without the
#     GE wage rise, which lowers every firm's value under free entry)
# and better_off/<id>_<grid>_heatmap.png: main5's reform heatmap (Functions/improve.jl plotPolHeatmaps), reform minus
# baseline by cash on hand and productivity, states with baseline mass above 1e-6.
using Dates, CSV, DataFrames, Printf
const SEARCH_DIR = @__DIR__
const MODEL_DIR  = normpath(joinpath(SEARCH_DIR, ".."))
const MODEL_FILE = joinpath(MODEL_DIR, "Baseline_noise.jl")
const DATA_DIR   = normpath(joinpath(MODEL_DIR, "..", "..", "Data Replication"))
ENV["GKSwstype"] = "100"
using LinearAlgebra
BLAS.set_num_threads(1)
let code = read(MODEL_FILE, String)
    i = findfirst("########## " * "CALIBRATION", code)
    task_local_storage(() -> include_string(Main, code[1:prevind(code, first(i))], MODEL_FILE), :SOURCE_PATH, MODEL_FILE)
end
cd(MODEL_DIR)
include(joinpath(MODEL_DIR, "Functions", "Results.jl"))

id = ARGS[1]
grid_over = Dict{Symbol, Float64}(Symbol(split(a, "=")[1]) => parse(Float64, split(a, "=")[2]) for a in ARGS[2:end])
run_dir = joinpath(SEARCH_DIR, "runs", id)
diag = first(CSV.read(joinpath(run_dir, "DIAG.csv"), DataFrame; types = String))
dv(k) = (k in names(diag) && !ismissing(diag[k]) && diag[k] != "") ? parse(Float64, diag[k]) : nothing

# external settings as in Search.jl: the recorded values, then the run's set: overrides (DIAG columns x_<name>)
ext = Dict{Symbol, Float64}()
for k in ("beta", "delta", "phi_a", "phi_af", "kappa", "rho_e", "sigma_e", "nul_e", "DRS", "size_by_debt", "d2c_cut", "tauchen_lo")
    v = dv(k); v === nothing || (ext[Symbol(k)] = v)
end
for k in names(diag)
    startswith(k, "x_") && dv(k) !== nothing && (ext[Symbol(k[3:end])] = dv(k))
end
merge!(ext, grid_over)
vals = load_calibration(joinpath(run_dir, "CALIB.csv"))
reset_params!(); set_params!((; ext...)); set_params!(vals)
g = gridsize()
gtag = "$(g.x_size)_$(g.e_size)_$(g.k_size)_$(g.b_size)"
println("$id on grid $gtag; externals ", ext); flush(stdout)

tag = diag["targets"]
data = read_targets(joinpath(DATA_DIR, "calib_targets_$tag.csv"))
share = data_value(data, "Share of small firms")

ge_try(v, c_e) = (for b in [(0.9, 1.3), (0.8, 1.6), (0.7, 2.2)]
                      try return ge_eq(v, c_e; bracket = b) catch e; println("ge_eq failed on $b") end
                  end; error("no GE"))

t0 = time()
B0 = solve_eq()
Fcut = fit_Fcut(B0.SumPol, B0.Fmat, B0.f0, share; zeta_R = parameters().zeta_Rl)
set_params!((Fcut = Fcut,))
B = solve_eq()
R = ge_try((zeta_Rl = parameters().zeta_Rl / 2,), B.c_e)
P = solve_eq((zeta_Rl = parameters().zeta_Rl / 2,); wage = B.wage)   # reform at the baseline wage: the credit effect alone
println("solved in ", round((time() - t0) / 60, digits = 1), " min"); flush(stdout)

# small-firm averages of the reform (or baseline) policies S over states with weights w; small set = baseline policy b' < Fcut
function small_avg(S, w, small)
    n = size(S, 1) - 1
    b, k, q = S[1:n, 4], S[1:n, 3], S[1:n, 9]
    bor = small .& (b .!= 0)
    rate_ok = bor .& (q .>= 0.5)
    d2c_ok = small .& (k .!= 0) .& (q .>= 0.5) .& (b ./ max.(k, 1e-12) .<= parameters().d2c_cut)
    liq_ok = bor .& (S[1:n, 8] .!= 0)
    avg(x, m) = sum(w[m] .* x[m]) / sum(w[m])
    (rate = 100 * avg(1 ./ q .- 1, rate_ok), d2c = avg(b ./ max.(k, 1e-12), d2c_ok), CFrel = avg(S[1:n, 17], bor),
     CFsh = sum(w[bor] .* b[bor] .* S[1:n, 17][bor]) / sum(w[bor] .* b[bor]), liq = avg(S[1:n, 14], liq_ok),
     value = avg(S[1:n, 18], small), mass = sum(w[small]))
end
n = size(B.SumPol, 1) - 1
smallB = (B.SumPol[1:n, 4] .< Fcut) .& (B.SumPol[1:n, 4] .!= 0)
smallR = (R.SumPol[1:n, 4] .< Fcut) .& (R.SumPol[1:n, 4] .!= 0)
base = small_avg(B.SumPol, B.mu[1:n], smallB)
tot  = small_avg(R.SumPol, R.mu[1:n], smallR)          # small firms in the reform equilibrium
same = small_avg(R.SumPol, B.mu[1:n], smallB)          # reform policies at the baseline's small-firm states
pe   = small_avg(P.SumPol, B.mu[1:n], smallB)          # the same, reform at the baseline wage

agg(E, r) = E.ss[r, "values"]
row = Pair{String, Any}["id" => id, "grid" => gtag, "Fcut" => Fcut, "minutes" => round((time() - t0) / 60, digits = 1),
    "d2c" => agg(B, "avg_b2a"), "rate" => agg(B, "avg intrate"), "liq" => agg(B, "avg_liqprob"), "CFrel" => agg(B, "CF reliance"),
    "CFsh" => agg(B, "CF share"),
    "d_d2c" => data_value(data, "Debt to Collateral"), "d_rate" => data_value(data, "Interest rate"), "d_liq" => data_value(data, "Liquidation prob."),
    "d_CFrel" => data_value(data, "CF-reliance"), "d_CFsh" => data_value(data, "CF-based to total debt")]
for (k, r) in ("liq" => "liqprob", "rate" => "intrate", "d2c" => "b2a", "CFrel" => "CF reliance", "CFsh" => "CF share")
    push!(row, "s_" * k => B.sme[r, "SME Values"], "l_" * k => B.sme[r, "LE Values"])
    push!(row, "r_s_" * k => R.sme[r, "SME Values"], "r_l_" * k => R.sme[r, "LE Values"])
end
push!(row, "YL" => 100 * (agg(R, "YtoL") / agg(B, "YtoL") - 1), "mass" => 100 * (agg(R, "totalmass") / agg(B, "totalmass") - 1), "wage" => R.wage)
for (lab, s) in ("b" => base, "t" => tot, "w" => same, "p" => pe), k in keys(s)
    push!(row, "sm_$(lab)_$(k)" => s[k])
end
out = joinpath(SEARCH_DIR, "better_off"); mkpath(out)
CSV.write(joinpath(out, "$(id)_$(gtag).csv"), DataFrame([k => [v] for (k, v) in row]))
@printf("%s %s: Y/L %+.2f%%; small rate %.2f -> %.2f (same states %.2f); d2c %.3f -> %.3f (same %.3f); V %.1f -> %.1f (same %.1f)\n",
        id, gtag, row[findfirst(p -> p.first == "YL", row)].second, base.rate, tot.rate, same.rate, base.d2c, tot.d2c, same.d2c,
        base.value, tot.value, same.value)

# main5's reform heatmap: reform minus baseline at each state (x, e), states with baseline mass above 1e-6
S, S0, mu = B.SumPol[1:n, :], R.SumPol[1:n, :], B.mu[1:n]
keep = (S[:, 1] .>= 0) .& (S[:, 1] .<= 50000) .& (S[:, 2] .>= 3) .& (mu .> 1e-6)
xs = sort(unique(S[(S[:, 1] .>= 0) .& (S[:, 1] .<= 50000), 1])); es = sort(unique(S[S[:, 2] .>= 3, 2]))
function hmat(v)
    M = fill(NaN, length(es), length(xs))
    for i in findall(keep)
        M[searchsortedfirst(es, S[i, 2]), searchsortedfirst(xs, S[i, 1])] = v[i]
    end
    M
end
b2a = S0[:, 4] ./ (S0[:, 3] .+ S0[:, 1]) .- S[:, 4] ./ (S[:, 3] .+ S[:, 1]); b2a[abs.(b2a) .> 0.3] .= NaN
dq = 100 .* (1 ./ S0[:, 9] .- 1 ./ S[:, 9]); dq[abs.(dq) .> 4] .= NaN
small_mark = hmat(Float64.(smallB))
hm(v, t) = heatmap(string.(round.(Int, xs)), string.(round.(es, digits = 1)), hmat(v), color = :viridis, xlabel = "Cash on Hand",
                   ylabel = "Productivity", title = t, xrotation = 45, tickfontsize = 7, titlefontsize = 11)
p = plot(hm(S0[:, 14] .- S[:, 14], "A: Liquidation Probability"), hm(S0[:, 17] .- S[:, 17], "B: CFL reliance"),
         hm(b2a, "C: Debt to Collateral"), hm(dq, "D: Interest Rate"),
         heatmap(string.(round.(Int, xs)), string.(round.(es, digits = 1)), small_mark, color = cgrad([:lightgray, :darkorange]),
                 colorbar = false, xlabel = "Cash on Hand", ylabel = "Productivity", title = "Small firms (baseline b' < Fcut)",
                 xrotation = 45, tickfontsize = 7, titlefontsize = 11),
         layout = (3, 2), size = (1000, 1050), margin = 4mm, plot_title = "$id, grid $gtag: reform minus baseline")
savefig(p, joinpath(out, "$(id)_$(gtag)_heatmap.png"))
println("saved ", joinpath(out, "$(id)_$(gtag).csv"))

# rate_decomp.jl <run id> [name=value ...]: where does the change in the small-firm interest rate in the zeta/2 reform come from?
# Same setup as better_off.jl (the run's CALIB.csv and externals, optional grid overrides). Splits the move from the
# baseline small-firm rate to the reform one (all small firms in each equilibrium) into
#   credit (reform policies at the baseline wage, baseline small states and weights), wage (the same at the GE wage),
#   reweighting (baseline small states, reform weights), and set change (states that become small or stop being small),
# and splits small firms by age: the stationary distribution is mu = m * sum_t M^t f0' (stat_dist), so cohort t is
# m * M^t f0'. Writes better_off/<id>_<grid>_ratedecomp.csv.
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

n = size(B.SumPol, 1) - 1
small(E) = (E.SumPol[1:n, 4] .< Fcut) .& (E.SumPol[1:n, 4] .!= 0)
sB, sR = small(B), small(R)
# the small-firm rate as in better_off.jl and sumSSsme: mean of 1/q - 1 over borrowers with q >= 0.5, weights w
function rate(S, w, s)
    q = S[1:n, 9]; ok = s .& (S[1:n, 4] .!= 0) .& (q .>= 0.5)
    (rate = 100 * sum(w[ok] .* (1 ./ q[ok] .- 1)) / sum(w[ok]), mass = sum(w[ok]))
end
wB, wR = B.mu[1:n], R.mu[1:n]
steps = [("baseline", rate(B.SumPol, wB, sB)),
         ("credit: reform at baseline wage, baseline small states and weights", rate(P.SumPol, wB, sB)),
         ("+ GE wage", rate(R.SumPol, wB, sB)),
         ("+ reform weights on baseline small states", rate(R.SumPol, wR, sB)),
         ("+ small set of the reform = reform small firms", rate(R.SumPol, wR, sR))]
newS, gone, both = sR .& .!sB, sB .& .!sR, sR .& sB
groups = [("reform small, small in both", rate(R.SumPol, wR, both)), ("reform small, newly small", rate(R.SumPol, wR, newS)),
          ("baseline small that stop being small (reform policies, reform weights)", rate(R.SumPol, wR, gone)),
          ("baseline small, small in both (baseline)", rate(B.SumPol, wB, both)), ("baseline small, stop being small (baseline)", rate(B.SumPol, wB, gone))]

# cohorts: mu = m * sum_t M^t f0', f0' = (I - diag(xpol)) f0, M = Fmat (I - diag(xpol))
function cohorts(E; T = 5)
    X = Diagonal(1 .- E.xpol); c = E.m .* (X * E.f0); M = E.Fmat * X
    out = Vector{Vector{Float64}}()
    for t in 0:T-1
        push!(out, c[1:n]); c = M * c
    end
    (age0 = out[1], age1_4 = sum(out[2:T]), age5p = E.mu[1:n] .- sum(out))
end
cB, cR = cohorts(B), cohorts(R)
ages = Tuple{String, Any}[]
for (lab, k) in (("age 0 (entrants)", :age0), ("age 1-4", :age1_4), ("age 5+", :age5p))
    push!(ages, ("baseline small, " * lab, rate(B.SumPol, getfield(cB, k), sB)))
    push!(ages, ("reform small, " * lab, rate(R.SumPol, getfield(cR, k), sR)))
    push!(ages, ("baseline small states with reform policies and baseline weights, " * lab, rate(R.SumPol, getfield(cB, k), sB)))
end
rows = DataFrame(part = String[], label = String[], rate = Float64[], mass = Float64[])
for (part, L) in (("steps", steps), ("groups", groups), ("ages", ages)), (l, r) in L
    push!(rows, (part, l, r.rate, r.mass))
end
println("check: cohorts add up to mu: ", maximum(abs.(cB.age0 .+ cB.age1_4 .+ cB.age5p .- wB)), " (age5p min ", minimum(cB.age5p), ")")
CSV.write(joinpath(SEARCH_DIR, "better_off", "$(id)_$(gtag)_ratedecomp.csv"), rows)
for r in eachrow(rows)
    @printf("%-7s %-95s %6.2f  mass %8.4f\n", r.part, r.label, r.rate, r.mass)
end

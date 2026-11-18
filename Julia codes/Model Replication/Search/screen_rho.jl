# screen_rho.jl <chunkfile>: screens the internal parameters under alternative productivity processes.
# Each line of <chunkfile> is a chunk "rho s j1 j2": persistence rho, stationary sd s of ln e, Halton points j1..j2 of
# the box below. sigma_e = s * sqrt(1 - rho^2); nul_e is re-centred so the mean of e^(1/(1-DRS)) stays at its value under
# parameters() (nul_e from screens/nul_*.csv, written by calib_nul.jl, when the cell is there: the entry value at the
# parameters() values of the internal parameters equals its baseline value). Every point records the aggregate fit (the five targeted moments, which are the same for every size cut)
# and the small and large firm moments under the $100m (L8) and $1bn (L9) ltq splits.
# Writes screens/rho<rho>_s<s>_<j1>.csv per chunk.
using Dates, CSV, DataFrames
const CHUNKFILE = abspath(ARGS[1])
const SEARCH_DIR = @__DIR__
const MODEL_DIR  = normpath(joinpath(SEARCH_DIR, ".."))
const MODEL_FILE = joinpath(MODEL_DIR, "Baseline.jl")
const DATA_DIR   = normpath(joinpath(MODEL_DIR, "..", "..", "Data Replication"))
ENV["GKSwstype"] = "100"
using LinearAlgebra
BLAS.set_num_threads(1)
let code = read(MODEL_FILE, String)
    i = findfirst("########## " * "CALIBRATION", code)
    task_local_storage(() -> include_string(Main, code[1:prevind(code, first(i))], MODEL_FILE), :SOURCE_PATH, MODEL_FILE)
end
cd(MODEL_DIR)

reset_params!()
const P0 = parameters()
const S0 = P0.sigma_e / sqrt(1 - P0.rho_e^2)
const MU0 = P0.nul_e
const splits = [("L8", read_targets(joinpath(DATA_DIR, "calib_targets_L8.csv"))), ("L9", read_targets(joinpath(DATA_DIR, "calib_targets_L9.csv")))]
const calib_moments = [("Debt to Collateral", "avg_b2a", "b2a", true), ("Interest rate", "avg intrate", "intrate", true),
                       ("Liquidation prob.", "avg_liqprob", "liqprob", true), ("CF-reliance", "CF reliance", "CF reliance", true),
                       ("CF-based to total debt", "CF share", "CF share", true)]
const box = (pc = (15.0, 60.0), pdef_exo_l = (0.003, 0.05), zeta_Rl = (50.0, 2500.0), phi_c_hh = (0.2, 0.95), phi_c = (0.05, 0.5))
const names_ = keys(box)
mkpath(joinpath(SEARCH_DIR, "screens"))
const NUL = Dict{Tuple{Float64, Float64}, Float64}()
for f in filter(x -> startswith(x, "nul_") && endswith(x, ".csv"), readdir(joinpath(SEARCH_DIR, "screens")))
    for r in eachrow(CSV.read(joinpath(SEARCH_DIR, "screens", f), DataFrame))
        NUL[(r.rho, r.s)] = r.nul_e
    end
end

for line in readlines(CHUNKFILE)
    isempty(strip(line)) && continue
    f = split(strip(line))
    rho = parse(Float64, f[1]); s = parse(Float64, f[2]); j1 = parse(Int, f[3]); j2 = parse(Int, f[4])
    ext = (rho_e = rho, sigma_e = s * sqrt(1 - rho^2), nul_e = get(NUL, (rho, s), MU0 + (S0^2 - s^2) / (2 * (1 - P0.DRS))))
    println(now(), " chunk rho $rho s $s points $j1-$j2: ", ext); flush(stdout)
    reset_params!(); set_params!(ext)
    rows = NamedTuple[]
    for j in j1:j2
        vals = NamedTuple{names_}(Tuple(box[k][1] + halton(j, HALTON_BASES[d]) * (box[k][2] - box[k][1]) for (d, k) in enumerate(names_)))
        vals.phi_c <= vals.phi_c_hh || continue          # valid: (1 - phi_c_hh) + phi_c <= 1
        t0 = time()
        r = try
            with_params(vals) do
                SumPol, e_chain, Fmat = FirmOptim(1.0)
                c_e, f0 = EntryValue(SumPol, e_chain)
                ss = sumSS(SumPol, Fmat, f0; zeta_R = parameters().zeta_Rl, rate_digits = nothing)
                t = targeted_table(ss, calib_moments, splits[1][2])
                out = (rho = rho, s = s, sigma_e = ext.sigma_e, nul_e = ext.nul_e, j = j)
                out = merge(out, vals, (error = sum(t.error), d2c = t.model[1], rate = t.model[2], liq = t.model[3], CFrel = t.model[4], CFsh = t.model[5]))
                for (tag, tg) in splits
                    Fcut = fit_Fcut(SumPol, Fmat, f0, data_value(tg, "Share of small firms"); zeta_R = parameters().zeta_Rl)
                    sme = with_params(() -> sumSSsme(SumPol, Fmat, f0; zeta_R = parameters().zeta_Rl, rate_digits = nothing), (Fcut = Fcut,))
                    g(row, col) = sme[row, col]
                    out = merge(out, NamedTuple{Tuple(Symbol.(["$(tag)_$k" for k in ("s_liq", "l_liq", "s_rate", "l_rate", "s_d2c", "l_d2c", "s_CFrel", "l_CFrel", "share")]))}(
                        (g("liqprob", "SME Values"), g("liqprob", "LE Values"), g("intrate", "SME Values"), g("intrate", "LE Values"),
                         g("b2a", "SME Values"), g("b2a", "LE Values"), g("CF reliance", "SME Values"), g("CF reliance", "LE Values"), g("share", "SME Values"))))
                end
                out
            end
        catch e
            println("point $j failed: ", sprint(showerror, e)); nothing
        end
        r === nothing && continue
        push!(rows, merge(r, (secs = round(time() - t0, digits = 1),)))
        println(now(), " point $j: error ", round(r.error, digits = 3), " liq ", round(r.liq, digits = 2), " L8 s/l ",
                round(r.L8_s_liq, digits = 2), "/", round(r.L8_l_liq, digits = 2), " (", round(time() - t0, digits = 1), " s)"); flush(stdout)
    end
    isempty(rows) || CSV.write(joinpath(SEARCH_DIR, "screens", "rho$(rho)_s$(s)_$(j1).csv"), DataFrame(rows))
end
println(now(), " done")

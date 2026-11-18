# l9_moments.jl <out.csv> <id>...: for finished runs (runs/<id>/CALIB.csv), the small and large firm moments under the
# $100m (L8) and $1bn (L9) ltq splits at the calibrated values, without the results stage. Externals rho_e, sigma_e,
# nul_e come from runs/<id>/DIAG.csv.
const OUT = abspath(ARGS[1])
const SEARCH_DIR = @__DIR__
const MODEL_FILE = normpath(joinpath(SEARCH_DIR, "..", "Baseline.jl"))
const DATA_DIR = normpath(joinpath(SEARCH_DIR, "..", "..", "..", "Data Replication"))
ENV["GKSwstype"] = "100"
using LinearAlgebra, CSV, DataFrames
BLAS.set_num_threads(1)
let code = read(MODEL_FILE, String)
    i = findfirst("########## " * "CALIBRATION", code)
    task_local_storage(() -> include_string(Main, code[1:prevind(code, first(i))], MODEL_FILE), :SOURCE_PATH, MODEL_FILE)
end
cd(dirname(MODEL_FILE))
splits = [("L8", read_targets(joinpath(DATA_DIR, "calib_targets_L8.csv"))), ("L9", read_targets(joinpath(DATA_DIR, "calib_targets_L9.csv")))]
rows = NamedTuple[]
for id in ARGS[2:end]
    c = CSV.read(joinpath(SEARCH_DIR, "runs", id, "CALIB.csv"), DataFrame; types = String)
    d = CSV.read(joinpath(SEARCH_DIR, "runs", id, "DIAG.csv"), DataFrame; types = String)
    iv(n) = parse(Float64, c[(c.section .== "internal parameter") .& (c.name .== n), :value][1])
    ext = (rho_e = parse(Float64, d.rho_e[1]), sigma_e = parse(Float64, d.sigma_e[1]), nul_e = parse(Float64, d.nul_e[1]), size_by_debt = 1.0)
    vals = (pc = iv("pc"), pdef_exo_l = iv("pdef_exo_l"), zeta_Rl = iv("zeta_Rl"), phi_c_hh = iv("phi_c_hh"), phi_c = iv("phi_c"))
    reset_params!(); set_params!(ext)
    out = with_params(vals) do
        SumPol, e_chain, Fmat = FirmOptim(1.0)
        c_e, f0 = EntryValue(SumPol, e_chain)
        o = (id = id, rho = ext.rho_e, s = ext.sigma_e / sqrt(1 - ext.rho_e^2))
        for (tag, tg) in splits
            Fcut = fit_Fcut(SumPol, Fmat, f0, data_value(tg, "Share of small firms"); zeta_R = parameters().zeta_Rl)
            sme = with_params(() -> sumSSsme(SumPol, Fmat, f0; zeta_R = parameters().zeta_Rl, rate_digits = nothing), (Fcut = Fcut,))
            g(r, col) = sme[r, col]
            ks = ("liqprob", "intrate", "b2a", "CF reliance", "CF share")
            o = merge(o, NamedTuple{Tuple(Symbol.(vcat(["$(tag)_s_$k" for k in ("liq", "rate", "d2c", "CFrel", "CFsh")], ["$(tag)_l_$k" for k in ("liq", "rate", "d2c", "CFrel", "CFsh")])))}(
                Tuple(vcat([g(k, "SME Values") for k in ks], [g(k, "LE Values") for k in ks]))))
        end
        o
    end
    push!(rows, out); println(out); flush(stdout)
    CSV.write(OUT, DataFrame(rows))
end

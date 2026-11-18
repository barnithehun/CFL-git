# sweep_param.jl <model file> <param> <out.csv> <value>...: the model file (in the model folder) at the parameters() values for each value of <param>:
# aggregate fit (targets of the $100m ltq split, L8) and small/large firm moments under the $100m (L8) and $1bn (L9) splits
const MODEL_NAME = ARGS[1]; const PNAME = Symbol(ARGS[2]); const OUT = abspath(ARGS[3])
const SEARCH_DIR = @__DIR__
const MODEL_FILE = normpath(joinpath(SEARCH_DIR, "..", MODEL_NAME))
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
cm = [("Debt to Collateral", "avg_b2a", "b2a", true), ("Interest rate", "avg intrate", "intrate", true),
      ("Liquidation prob.", "avg_liqprob", "liqprob", true), ("CF-reliance", "CF reliance", "CF reliance", true),
      ("CF-based to total debt", "CF share", "CF share", true)]
rows = NamedTuple[]
for sv in parse.(Float64, ARGS[4:end])
    reset_params!()
    r = with_params(NamedTuple{(PNAME,)}((sv,))) do
        SumPol, e_chain, Fmat = FirmOptim(1.0)
        c_e, f0 = EntryValue(SumPol, e_chain)
        ss = sumSS(SumPol, Fmat, f0; zeta_R = parameters().zeta_Rl, rate_digits = nothing)
        t = targeted_table(ss, cm, splits[1][2])
        o = (value = sv, colsum_err = maximum(abs.(sum(Fmat, dims = 1)[1:end-1] .- 1)), c_e = c_e, error = sum(t.error), d2c = t.model[1], rate = t.model[2], liq = t.model[3], CFrel = t.model[4], CFsh = t.model[5])
        for (tag, tg) in splits
            Fcut = fit_Fcut(SumPol, Fmat, f0, data_value(tg, "Share of small firms"); zeta_R = parameters().zeta_Rl)
            sme = with_params(() -> sumSSsme(SumPol, Fmat, f0; zeta_R = parameters().zeta_Rl, rate_digits = nothing), (Fcut = Fcut,))
            o = merge(o, NamedTuple{Tuple(Symbol.(["$(tag)_s_liq", "$(tag)_l_liq", "$(tag)_s_rate", "$(tag)_l_rate", "$(tag)_s_CFrel", "$(tag)_l_CFrel"]))}(
                (sme["liqprob", "SME Values"], sme["liqprob", "LE Values"], sme["intrate", "SME Values"], sme["intrate", "LE Values"],
                 sme["CF reliance", "SME Values"], sme["CF reliance", "LE Values"])))
        end
        o
    end
    push!(rows, r)
    println("$PNAME $sv: Fmat colsum err ", round(r.colsum_err, sigdigits = 2), " error ", round(r.error, digits = 3), " liq ", round(r.liq, digits = 3), " L8 s/l ", round(r.L8_s_liq, digits = 2), "/", round(r.L8_l_liq, digits = 2),
            " L9 s/l ", round(r.L9_s_liq, digits = 2), "/", round(r.L9_l_liq, digits = 2), " rate ", round(r.rate, digits = 2), " CFsh ", round(r.CFsh, digits = 2), " c_e ", round(r.c_e, digits = 1)); flush(stdout)
    CSV.write(OUT, DataFrame(rows))
end

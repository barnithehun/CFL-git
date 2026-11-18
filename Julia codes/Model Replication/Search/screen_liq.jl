# screen_liq.jl <tag> <n> [name=value ...]: the aggregate fit and the small and large firm moments at n Halton points of a
# wide box of the internal parameters, on the targets file of <tag> (as in Search.jl). Answers where in the parameter space
# large firms are liquidated at all. Writes screen_<tag>_<n>.csv in this folder. name=value sets external overrides.
using Dates, CSV, DataFrames
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

tag = ARGS[1]; n = parse(Int, ARGS[2])
# name=value: external override; name=lo:hi: screening range of an internal parameter
ext = Dict{Symbol, Float64}(Symbol(split(a, "=")[1]) => parse(Float64, split(a, "=")[2]) for a in ARGS[3:end] if !occursin(":", a))
box_over = Dict{Symbol, Tuple{Float64, Float64}}(Symbol(split(a, "=")[1]) => Tuple(parse.(Float64, split(split(a, "=")[2], ":"))) for a in ARGS[3:end] if occursin(":", a))
haskey(ext, :phi_a) && (ext[:phi_af] = ext[:phi_a])
tfile = joinpath(DATA_DIR, tag == "T0" ? "calib_targets.csv" : "calib_targets_$tag.csv")
data_targets = read_targets(tfile)
share = data_value(data_targets, "Share of small firms")
calib_moments = [("Debt to Collateral", "avg_b2a", "b2a", true), ("Interest rate", "avg intrate", "intrate", true),
                 ("Liquidation prob.", "avg_liqprob", "liqprob", true), ("CF-reliance", "CF reliance", "CF reliance", true),
                 ("CF-based to total debt", "CF share", "CF share", true)]
box = (pc = (15.0, 60.0), pdef_exo_l = (0.003, 0.05), zeta_Rl = (50.0, 1500.0), phi_c_hh = (0.2, 0.9), phi_c = (0.05, 0.5))
box = merge(box, (; box_over...))
names_ = keys(box)

rows = NamedTuple[]
reset_params!(); set_params!((; ext...))
for j in 1:n
    vals = NamedTuple{names_}(Tuple(box[k][1] + halton(j, HALTON_BASES[d]) * (box[k][2] - box[k][1]) for (d, k) in enumerate(names_)))
    vals.phi_c <= vals.phi_c_hh || continue          # valid: (1 - phi_c_hh) + phi_c <= 1
    r = try
        with_params(vals) do
            SumPol, e_chain, Fmat = FirmOptim(1.0)
            c_e, f0 = EntryValue(SumPol, e_chain)
            ss = sumSS(SumPol, Fmat, f0; zeta_R = parameters().zeta_Rl, rate_digits = nothing)
            t = targeted_table(ss, calib_moments, data_targets)
            Fcut = fit_Fcut(SumPol, Fmat, f0, share; zeta_R = parameters().zeta_Rl)
            sme = with_params(() -> sumSSsme(SumPol, Fmat, f0; zeta_R = parameters().zeta_Rl, rate_digits = nothing), (Fcut = Fcut,))
            merge(vals, (error = sum(t.error), d2c = t.model[1], rate = t.model[2], liq = t.model[3], CFrel = t.model[4], CFsh = t.model[5],
                         s_liq = sme["liqprob", "SME Values"], l_liq = sme["liqprob", "LE Values"],
                         s_rate = sme["intrate", "SME Values"], l_rate = sme["intrate", "LE Values"],
                         s_CFrel = sme["CF reliance", "SME Values"], l_CFrel = sme["CF reliance", "LE Values"],
                         s_d2c = sme["b2a", "SME Values"], l_d2c = sme["b2a", "LE Values"], share = sme["share", "SME Values"]))
        end
    catch e
        println("point $j failed: ", sprint(showerror, e)); nothing
    end
    r === nothing || push!(rows, r)
    j % 20 == 0 && (println("$j points"); flush(stdout))
end
out = joinpath(SEARCH_DIR, "screen_$(tag)_$(n)$(isempty(ext) ? "" : "_" * join(["$(k)$(v)" for (k, v) in ext if k != :phi_af], "_"))$(isempty(box_over) ? "" : "_box").csv")
CSV.write(out, sort(DataFrame(rows), :error))
println("saved ", out)

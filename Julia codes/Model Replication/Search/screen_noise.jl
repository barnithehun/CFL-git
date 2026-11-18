# screen_noise.jl <tags> <n> <worker> <workers> [name=value ...] [name=lo:hi ...]: the noise model (Baseline_noise.jl) at
# Halton points j = worker+1, worker+1+workers, ... <= n of a wide box of the internal parameters (sigma_v included), on the
# targets file of each tag in <tags> (comma-separated; the first names the output). Records the five aggregate moments and the small and large firm moments at the data's share of
# small firms. name=value: external override (e.g. tau_tie=2, phi_a=0.425); name=lo:hi: the box of an internal parameter
# (sv_slope=-0.5:0.5 adds the size-dependent noise). The model is solved once per point; the size moments are computed
# for each tag (columns <moment>_<tag>). Writes screens/noise_<label>_<n>_w<worker>.csv, label = <first tag>[_sv].
using Dates, CSV, DataFrames
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

tags = split(ARGS[1], ","); n = parse(Int, ARGS[2]); worker = parse(Int, ARGS[3]); workers = parse(Int, ARGS[4])
rest = ARGS[5:end]
ext = Dict{Symbol, Float64}(Symbol(split(a, "=")[1]) => parse(Float64, split(a, "=")[2]) for a in rest if !occursin(":", a))
box_over = Dict{Symbol, Tuple{Float64, Float64}}(Symbol(split(a, "=")[1]) => Tuple(parse.(Float64, split(split(a, "=")[2], ":"))) for a in rest if occursin(":", a))
haskey(ext, :phi_a) && (ext[:phi_af] = ext[:phi_a])
tfile(tag) = joinpath(DATA_DIR, tag in ("", "T0") ? "calib_targets.csv" : "calib_targets_$tag.csv")
targets = Dict(tag => read_targets(tfile(tag)) for tag in tags)
data_targets = targets[tags[1]]
by_assets = Dict(tag => occursin("size_var,assets", read(replace(tfile(tag), ".csv" => "_settings.csv"), String)) for tag in tags)
println("assets splits: ", [t for t in tags if by_assets[t]])
calib_moments = [("Debt to Collateral", "avg_b2a", "b2a", true), ("Interest rate", "avg intrate", "intrate", true),
                 ("Liquidation prob.", "avg_liqprob", "liqprob", true), ("CF-reliance", "CF reliance", "CF reliance", true),
                 ("CF-based to total debt", "CF share", "CF share", true)]
box = (pc = (10.0, 120.0), pdef_exo_l = (0.001, 0.05), zeta_Rl = (100.0, 3000.0), phi_c_hh = (0.2, 1.0), phi_c = (0.05, 0.6),
       sigma_v = (0.3, 2.5))
box = merge(box, (; box_over...))
names_ = keys(box)

out = joinpath(SEARCH_DIR, "screens", "noise_$(tags[1])$(haskey(box_over, :sv_slope) ? "_sv" : "")_$(n)_w$(worker).csv")
rows = NamedTuple[]
reset_params!(); set_params!((; ext...))
t0 = time()
for j in (worker + 1):workers:n
    vals = NamedTuple{names_}(Tuple(box[k][1] + halton(j, HALTON_BASES[d]) * (box[k][2] - box[k][1]) for (d, k) in enumerate(names_)))
    (1 - vals.phi_c_hh) + vals.phi_c <= 1 || continue          # valid
    r = try
        with_params(vals) do
            SumPol, e_chain, Fmat = FirmOptim(1.0)
            c_e, f0 = EntryValue(SumPol, e_chain)
            ss = sumSS(SumPol, Fmat, f0; zeta_R = parameters().zeta_Rl, rate_digits = nothing)
            t = targeted_table(ss, calib_moments, data_targets)
            sz = NamedTuple[]
            for tag in tags
                share = data_value(targets[tag], "Share of small firms")
                sbd = (size_by_debt = by_assets[tag] ? 0.0 : 1.0,)        # assets split: small by k', as Search.jl
                Fcut = with_params(() -> fit_Fcut(SumPol, Fmat, f0, share; zeta_R = parameters().zeta_Rl), sbd)
                sme = with_params(() -> sumSSsme(SumPol, Fmat, f0; zeta_R = parameters().zeta_Rl, rate_digits = nothing), merge(sbd, (Fcut = Fcut,)))
                v = (s_liq = sme["liqprob", "SME Values"], l_liq = sme["liqprob", "LE Values"],
                     s_rate = sme["intrate", "SME Values"], l_rate = sme["intrate", "LE Values"],
                     s_CFrel = sme["CF reliance", "SME Values"], l_CFrel = sme["CF reliance", "LE Values"],
                     s_d2c = sme["b2a", "SME Values"], l_d2c = sme["b2a", "LE Values"], share = sme["share", "SME Values"], Fcut = Fcut)
                push!(sz, NamedTuple{Tuple(Symbol(k, "_", tag) for k in keys(v))}(values(v)))
            end
            merge((point = j,), vals, (error = sum(t.error), c_e = c_e, d2c = t.model[1], rate = t.model[2], liq = t.model[3],
                         CFrel = t.model[4], CFsh = t.model[5]), sz...)
        end
    catch e
        println("point $j failed: ", sprint(showerror, e)); nothing
    end
    r === nothing || push!(rows, r)
    if length(rows) % 5 == 0 && !isempty(rows)
        CSV.write(out, DataFrame(rows))
        println("$(length(rows)) points, ", round((time() - t0) / 60, digits = 1), " min"); flush(stdout)
    end
end
CSV.write(out, DataFrame(rows))
println("saved ", out, " (", length(rows), " points)")

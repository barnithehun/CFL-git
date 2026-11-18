############ E00: Model_ext.jl against Baseline_noise46.jl ############
# Run twice: plain (Model_ext.jl through common.jl) and with ENV CHECK_ORIGINAL = 1 (Baseline_noise46.jl itself).
# The printed lines must be identical: the baseline, and the zeta/3 reform at a fixed wage of 1.01.
if get(ENV, "CHECK_ORIGINAL", "0") == "1"
    const MODEL_DIR = dirname(@__DIR__)
    let file = joinpath(MODEL_DIR, "Baseline_noise46.jl"), code = read(file, String)
        i = findfirst("########## " * "CALIBRATION", code)
        task_local_storage(() -> include_string(Main, code[1:first(i)-1], file), :SOURCE_PATH, file)
    end
    include(joinpath(MODEL_DIR, "Functions", "Results.jl"))
else
    include(joinpath(@__DIR__, "common.jl"))
end

reset_params!()
for (label, vals, wage) in [("baseline", NamedTuple(), 1.0), ("zeta/3 at w = 1.01", (zeta_Rl = parameters().zeta_Rl / 3,), 1.01)]
    t = @elapsed E = solve_eq(vals; wage = wage)
    println("CHECK ", label, ": SumPol ", hash(E.SumPol), " Fmat ", hash(Matrix(E.Fmat)), " mu ", hash(vec(E.mu)),
            " c_e ", E.c_e, " Y/L ", E.ss["YtoL", "values"], " liq s/l ", E.sme["liqprob", "SME Values"], " ", E.sme["liqprob", "LE Values"])
    println("TIME ", label, ": ", round(t, digits = 1), " s")
end

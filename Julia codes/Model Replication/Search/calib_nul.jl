# calib_nul.jl <out.csv> <rho> <s> [<rho> <s> ...]: for each productivity process (persistence rho, stationary sd s of ln e,
# sigma_e = s * sqrt(1 - rho^2)), finds the mean nul_e of ln e at which the entry value at the parameters() values of
# the internal parameters equals its value under parameters() (c_e0), by bisection. This keeps the scale of the economy
# (the entry cost) fixed when the process changes. Writes rho, s, sigma_e, nul_e, c_e to <out.csv>.
const OUT = abspath(ARGS[1])
const MODEL_FILE = normpath(joinpath(@__DIR__, "..", "Baseline.jl"))
ENV["GKSwstype"] = "100"
using LinearAlgebra, CSV, DataFrames
BLAS.set_num_threads(1)
let code = read(MODEL_FILE, String)
    i = findfirst("########## " * "CALIBRATION", code)
    task_local_storage(() -> include_string(Main, code[1:prevind(code, first(i))], MODEL_FILE), :SOURCE_PATH, MODEL_FILE)
end
cd(dirname(MODEL_FILE))
reset_params!()
const P0 = parameters()
ce(ext) = (reset_params!(); set_params!(ext); S = FirmOptim(1.0); EntryValue(S[1], S[2])[1])
const CE0 = ce((rho_e = P0.rho_e,))
println("c_e0 = ", CE0); flush(stdout)
rows = NamedTuple[]
for i in 2:2:length(ARGS)
    rho = parse(Float64, ARGS[i]); s = parse(Float64, ARGS[i + 1]); sig = s * sqrt(1 - rho^2)
    lo, hi = P0.nul_e - 1.0, P0.nul_e + 3.0
    f(m) = ce((rho_e = rho, sigma_e = sig, nul_e = m)) - CE0
    flo, fhi = f(lo), f(hi)
    (flo < 0 < fhi) || println("rho $rho s $s: no bracket, f(lo) = $flo, f(hi) = $fhi")
    for it in 1:16
        mid = (lo + hi) / 2; fm = f(mid)
        fm > 0 ? (hi = mid) : (lo = mid)
    end
    m = (lo + hi) / 2; c = f(m) + CE0
    println("rho $rho s $s sigma $sig: nul_e = $m, c_e = $c"); flush(stdout)
    push!(rows, (rho = rho, s = s, sigma_e = sig, nul_e = m, c_e = c))
    CSV.write(OUT, DataFrame(rows))
end

############ H1 (heavy, isolated): E16 with a refit at each phi_a ############
# E16 changes phi_a and keeps every other parameter. Here the six internal parameters are refitted (calibrate_lazy,
# the targets and settings of Baseline_noise46.jl's calibration block) at each phi_a, so every economy matches the
# same data moments, and then the zeta/3 reform is solved in GE. One lazy fit takes roughly half an hour to an hour.
# ENV H1_PHIS sets the values (comma separated; default 0.3,0.55), so each value can run in its own process.
include(joinpath(@__DIR__, "common.jl"))
dir = outdir("H1")

data_targets = read_targets(joinpath(MODEL_DIR, "..", "..", "Data Replication", "calib_targets_L9.csv"))
calib_moments = [
    ("Debt to Collateral",     "avg_b2a",     "b2a",         true),
    ("Interest rate",          "avg intrate", "intrate",     true),
    ("Liquidation prob.",      "avg_liqprob", "liqprob",     true),
    ("CF-reliance",            "CF reliance", "CF reliance", true),
    ("CF-based to total debt", "CF share",    "CF share",    true),
    ("Liquidation prob.",      "gap",         "liqprob",     true),
    ("CF-reliance",            "gap",         "CF reliance", true),
]
calib_params = (
    pc         = (parameters().pc,         0.0, Inf),
    pdef_exo_l = (parameters().pdef_exo_l, 0.0, 1.0),
    zeta_Rl    = (parameters().zeta_Rl,    0.0, Inf),
    phi_c_hh   = (parameters().phi_c_hh,   0.0, 1.0),
    phi_c      = (parameters().phi_c,      0.0, 1.0),
    sigma_v    = (parameters().sigma_v,    0.0, 5.0),
)
valid = p -> (1 - p.phi_c_hh) + p.phi_c <= 1
options = Optim.Options(iterations = 200, time_limit = 3 * 3600.0, f_reltol = 0.002, show_trace = false, store_trace = true)

phis = parse.(Float64, split(get(ENV, "H1_PHIS", "0.3,0.55"), ","))
for phi in phis
    fixed = (phi_a = phi, phi_af = phi)
    t = @elapsed Cal = with_params(() -> calibrate_lazy(calib_params, calib_moments, data_targets; options = options, valid = valid,
                                                        file = joinpath(dir, "CALIB_phi_a_$(phi).csv")), fixed)
    vals = merge(Cal.values, fixed)
    Bp = solve_eq(vals)                                              # wage 1, its own entry cost
    Rp = ge_eq(merge(vals, (zeta_Rl = vals.zeta_Rl / 3,)), Bp.c_e; bracket = (0.9, 1.3))
    ob, o = with_params(() -> outcomes(Bp), vals), with_params(() -> outcomes(Rp), vals)
    save_table(dir, "refit_phi_a_$(phi).csv", ["value"],
               ["fit minutes" => [t / 60], "Y/L gain of zeta/3, %" => [pct(Rp.wage, Bp.wage)], "C gain, %" => [pct(o["C"], ob["C"])],
                "C gain, main5 entry, %" => [pct(o["C main5"], ob["C main5"])], "baseline CF reliance" => [ob["CF reliance"]],
                "baseline liq share" => [100 * ob["liq share"]], "zeta" => [vals.zeta_Rl], "theta" => [vals.phi_c], "sigma_v" => [vals.sigma_v],
                "pc" => [vals.pc], "P_x" => [vals.pdef_exo_l], "phi_c_hh" => [vals.phi_c_hh]])
end

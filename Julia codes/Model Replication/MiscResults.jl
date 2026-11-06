###########################################################################
############ MISCELLANEOUS RESULTS: not in the paper, kept  ###############
###########################################################################
# Diagnostics and experiments from FinalVers/ that the paper does not use. Each block names the original code.
#
# How to run:
#  - Once, from baseline_calib.jl: run_misc = true at the end of its RESULTS section includes this file. It uses the
#    baseline and the reform solved there (BaseEq, ReformEq).
#  - On its own: first run baseline_calib.jl up to the CALIBRATION block (the model). The SETUP block below then solves
#    the baseline and the reform again from the newest results_outputs/<date>/RESULTS.csv (a few seconds), or from
#    results_file if you set it.
#  - Block by block: after SETUP, each block runs on its own.
# Figures are saved to results_outputs/misc_<date>/.

cd(@__DIR__)

########## SETUP ##########
@isdefined(FirmOptim) || error("run baseline_calib.jl up to the CALIBRATION block first: MiscResults.jl needs the model")
@isdefined(solve_eq) || include("Functions/Results.jl")
@isdefined(CapAlloc2) || include("Functions/CapAlloc2.jl")

if !(@isdefined(BaseEq) && @isdefined(ReformEq))
    results_file = @isdefined(results_file) ? results_file : newest_results()
    println("Solving the baseline and the reform saved in ", results_file)
    saved = load_eqs(results_file)
    reset_params!(); set_params!(saved["Baseline"].values)
    BaseEq   = solve_eq(; wage = saved["Baseline"].wage)
    ReformEq = solve_eq(saved["Reform"].values; wage = saved["Reform"].wage)
end
reset_params!(); set_params!(BaseEq.values)     # parameters() returns the baseline's values, which the functions below read
misc_dir = joinpath("results_outputs", "misc_" * Dates.format(now(), "yyyy-mm-dd_HHMM"))
mkpath(misc_dir)
savemisc(fig, name) = (savefig(fig, joinpath(misc_dir, "$name.png")); fig)

# The original functions read these globals: the baseline without a suffix, the reform with suffix 0
SumPol, e_chain, Fmat, f0, mu = BaseEq.SumPol, BaseEq.e_chain, BaseEq.Fmat, BaseEq.f0, BaseEq.mu
SumPol0, e_chain0, Fmat0, f00, mu0 = ReformEq.SumPol, ReformEq.e_chain, ReformEq.Fmat, ReformEq.f0, ReformEq.mu
wage, phi_c, zeta_Rl = BaseEq.wage, parameters().phi_c, parameters().zeta_Rl
_, e_size, _, _ = gridsize()
c_e_baseline = BaseEq.c_e


########## M1. Untargeted moments across firm sizes (Baseline.jl lines 571-580: PolperSize.jl, CapAlloc.jl) ##########
binnum = 10
# Pdef and CFL reliance; in calculating tau, this contains firms with Pdef = 0
savemisc(Ushape(binnum, SumPol, mu), "M1_Ushape")
# Gamma and CFL reliance; in calculating tau, this does not contain firms with Pdef = 0
savemisc(Xcross(binnum, SumPol, mu), "M1_Xcross")
# Average Debt policy and Interest rate across firm sizes
savemisc(QBplot(binnum, SumPol, SumPol0, mu, mu0), "M1_QBplot")
# Capital Allocation
savemisc(CapAlloc(SumPol, SumPol0, mu, mu0), "M1_CapAlloc")


########## M2. Debt schedule against leverage (Baseline.jl lines 583-588: DebtScedule in GamPol.jl) ##########
# Not run: DebtScedule and GridMake.jl follow an older model version (the old liquidation rule, the lender paying
# zeta in Pi_reo) and unpack parameters() in an old order, so their output does not describe this model.
# p1, p2 = DebtScedule(SumPol, 25, 17; phi_c = phi_c)
# p3, p4 = DebtScedule(SumPol, 25, 20; phi_c = phi_c)
# p5, p6 = DebtScedule(SumPol, 25, 22; phi_c = phi_c)
# plot(p3, p4, p5, p6, layout = (2, 2), size = (900, 700))


########## M3. Lifecycle simulations and policies (Baseline.jl lines 590-599: dynsim2.jl, plotPol.jl) ##########
# shows the average policies of firms of a given productivity across their lifecycle (random draws)
savemisc(dynsim2(e_size - 6, 10000), "M3_dynsim2_e$(e_size - 6)")
savemisc(dynsim2(e_size - 3, 10000), "M3_dynsim2_e$(e_size - 3)")
savemisc(dynsim2(e_size - 1, 10000), "M3_dynsim2_e$(e_size - 1)")

# shows optimal policies of firms with a certain productivity, across x-sizes
savemisc(plotPol(SumPol0, SumPol, e_size-6), "M3_plotPol")


########## M4. Firm states and policies in the stationary equilibrium (Baseline.jl lines 601-609: StatDist_plot.jl) ##########
binnum = 20
savemisc(plot(plotPDF(binnum, 'k', SumPol), plotPDF(binnum, 'b', SumPol), plotPDF(binnum, 'l', SumPol),
    plotPDF(binnum, 'y', SumPol), plotPDF(binnum, 'p', SumPol), plotPDF(binnum, 'v', SumPol), layout=(2,3), size=(1200, 800)), "M4_PDF")
savemisc(plot(plotCDF(binnum, 'k', SumPol), plotCDF(binnum, 'b', SumPol), plotCDF(binnum, 'l', SumPol),
    plotCDF(binnum, 'y', SumPol), plotCDF(binnum, 'p', SumPol), plotCDF(binnum, 'v', SumPol), layout=(2,3), size=(1200, 800)), "M4_CDF")

savemisc(plotXE(SumPol, mu, e_chain), "M4_XE")
savemisc(plotTauDist(SumPol), "M4_TauDist")
# production across productivity states in the reform (StatDist_plot.jl, not called in Baseline.jl)
savemisc(ProdDist(SumPol, mu, SumPol0, mu0), "M4_ProdDist")
# PrintPolOld(SumPol, mu)    # writes every state's policies to CFL_<date>.xlsx in this folder (ABcase.jl)


########## M5. Reform effects across states (appendix reform.png; Baseline.jl lines 616-620: improve.jl) ##########
# Including improve.jl also runs its line plotPolHeatmaps(SumPol, SumPol0, mu, 1)
include("Functions/improve.jl")
reformplot = savemisc(plotPolHeatmaps(SumPol, SumPol0, mu, 1), "M5_reform")


########## M6. Capital allocation and firm distribution, baseline and reform (CapAlloc2.jl's plots) ##########
Decomp = CapAlloc2(BaseEq, ReformEq)
savemisc(Decomp.capital, "M6_capital_allocation")
savemisc(Decomp.bar, "M6_decomposition")
savemisc(Decomp.distribution, "M6_firm_distribution")
print_table("Decomposition with CapAlloc2.jl's -0.003 / +0.003 adjustment of the extensive and policy margins",
            CapAlloc2(BaseEq, ReformEq; adjust = 0.003).table)


########## M7. Optimal policies over different reorganization costs (Baseline.jl lines 623-626: ZetaGam.jl, ZetaGamPlot.jl) ##########
# ZetaGam.jl called a 4-argument FirmOptim at old parameter values and its own grid; this loop solves the model at the
# baseline's values for each zeta (for all firms), about 1 s each.
zeta_R_vec = 0:1000:8000
SumPolZeta = zeros(size(SumPol)..., length(zeta_R_vec))
for (ind, zeta_RR) in enumerate(zeta_R_vec)
    println("Iteration: $ind, Value: $zeta_RR")
    SumPolZeta[:, :, ind] = FirmOptim(wage, phi_c, zeta_RR)[1]
end
savemisc(ZetaGamPlot(SumPolZeta, zeta_R_vec, 23, 15), "M7_ZetaGam")


########## M8. Targeted vs universal cost cuts (SmallCosts.jl lines 545-653, smallcosts.png) ##########
# At the baseline wage, as in SmallCosts.jl. The original grid 200:200:2600 is scaled to the baseline zeta.
zeta_Rs_vec = round.(zeta_Rl .* (0.1:0.1:1.0), digits = 2)
SC = Dict(k => zeros(length(zeta_Rs_vec), 2) for k in (:prodimp, :liqprob, :d2c, :CFrel, :intrate, :liqprob_SME, :liqprob_LE))
for (ind_s, zeta_Rs) in enumerate(zeta_Rs_vec)
    # targeted policy (column 1) and universal policy (column 2)
    for (col, vals) in enumerate(((zeta_Rs = zeta_Rs,), (zeta_Rl = zeta_Rs,)))
        eq = solve_eq(vals; wage = wage)
        SC[:prodimp][ind_s, col] = eq.c_e / c_e_baseline
        SC[:liqprob][ind_s, col] = eq.ss["avg_liqprob", "values"]
        SC[:d2c][ind_s, col]     = eq.ss["avg_b2a", "values"]
        SC[:CFrel][ind_s, col]   = eq.ss["CF reliance", "values"]
        SC[:intrate][ind_s, col] = eq.ss["avg intrate", "values"]
        SC[:liqprob_SME][ind_s, col] = eq.sme["liqprob", "SME Values"]
        SC[:liqprob_LE][ind_s, col]  = eq.sme["liqprob", "LE Values"]
    end
end

using StatsPlots
colors = [:lightsalmon :deepskyblue2]  # First column red, second column blue
label = ["Targeted" "Universal"]
xticks_labels = string.(round.(Int, zeta_Rs_vec))
# Function to get y-limits dynamically
function get_ylims(data)
    ymin, ymax = extrema(data)  # Get min and max of the data
    return (ymin - 0.1 * abs(ymax - ymin), ymax + 0.1 * abs(ymax - ymin))  # Add some padding
end
scbar(data, title) = groupedbar(data, bar_position=:dodge, title=title, labels=label, xticks=(1:length(zeta_Rs_vec), xticks_labels),
                                fill=colors, ylims=get_ylims(data), xrotation=45, xlabelfontsize=12, ylabelfontsize=12, tickfontsize=9, titlefontsize=13)
finplot = plot(scbar(SC[:liqprob], "A: Liquidation Probability"), scbar(SC[:CFrel], "C: CF reliance"),
               scbar(SC[:intrate], "D: Interest Rate"), scbar(SC[:d2c], "B: Debt to Collateral"), layout=(2,2), size=(900,600))
savemisc(finplot, "M8_smallcosts")


########## M9. Productivity over small and large firms' costs (improve.jl lines 1-84, commented out there) ##########
# At the baseline wage; c_e relative to the baseline, and small (SME) and large (LE) firm moments, on a grid of zeta_Rl x zeta_Rs.
zeta_Rl_vec = round.(zeta_Rl .* (0:1/6:1), digits = 2)
HM = Dict(k => zeros(length(zeta_Rl_vec), length(zeta_Rl_vec)) for k in (:prodimp, :liqprob_SME, :liqprob_LE, :d2c_SME, :d2c_LE, :CFrel_SME, :CFrel_LE))
for (ind_l, zeta_l) in enumerate(zeta_Rl_vec), (ind_s, zeta_s) in enumerate(zeta_Rl_vec)
    eq = solve_eq((zeta_Rl = zeta_l, zeta_Rs = zeta_s); wage = wage)
    HM[:prodimp][ind_l, ind_s] = eq.c_e / c_e_baseline
    for (k, row) in ((:liqprob, "liqprob"), (:d2c, "b2a"), (:CFrel, "CF reliance"))
        HM[Symbol(k, :_SME)][ind_l, ind_s] = eq.sme[row, "SME Values"]
        HM[Symbol(k, :_LE)][ind_l, ind_s]  = eq.sme[row, "LE Values"]
    end
end
hm(k, title) = heatmap(zeta_Rl_vec, zeta_Rl_vec, HM[k], xlabel = "zeta_Rs (small firms)", ylabel = "zeta_Rl (large firms)", title = title)
savemisc(hm(:prodimp, "Entry value relative to the baseline"), "M9_productivity_heatmap")
savemisc(plot(hm(:liqprob_SME, "Liq. prob., small"), hm(:liqprob_LE, "Liq. prob., large"), hm(:d2c_SME, "Debt/collateral, small"),
              hm(:d2c_LE, "Debt/collateral, large"), hm(:CFrel_SME, "CF reliance, small"), hm(:CFrel_LE, "CF reliance, large"),
              layout = (3, 2), size = (1000, 1100)), "M9_combined_heatmaps")


########## M10. Shifting costs: zeta and theta both cut by half (Baseline.jl lines 688-711, ShiftingCosts.jl) ##########
# General equilibrium, as in Baseline.jl's block; theta cut for all firms, and for small firms only (ShiftingCosts.jl's target_flag)
ShiftAllEq   = ge_eq((zeta_Rl = zeta_Rl / 2, phi_c = phi_c / 2), c_e_baseline)
ShiftSmallEq = ge_eq((zeta_Rl = zeta_Rl / 2, phi_c_s = phi_c / 2), c_e_baseline)
print_table("Shifting costs: changes from the baseline", change_table(["zeta/2 only" => (BaseEq, ReformEq),
    "zeta/2, theta/2 all firms" => (BaseEq, ShiftAllEq), "zeta/2, theta/2 small firms" => (BaseEq, ShiftSmallEq)]))


########## M11. Lower theta at the reform's zeta (improve.jl lines 153-210, commented out there) ##########
# At the baseline wage and zeta/2: theta - shift for all firms
shift_vec = 0:0.02:0.2
SH = Dict(k => zeros(length(shift_vec)) for k in (:prodimp, :liqprob, :d2c, :CFrel, :intrate, :liqprob_SME, :liqprob_LE, :d2c_SME, :d2c_LE, :CFrel_SME, :CFrel_LE))
for (ind, shift) in enumerate(shift_vec)
    eq = solve_eq((zeta_Rl = zeta_Rl / 2, phi_c = phi_c - shift); wage = wage)
    SH[:prodimp][ind] = eq.c_e / c_e_baseline
    SH[:liqprob][ind] = eq.ss["avg_liqprob", "values"]
    SH[:d2c][ind]     = eq.ss["avg_b2a", "values"]
    SH[:CFrel][ind]   = eq.ss["CF reliance", "values"]
    SH[:intrate][ind] = eq.ss["avg intrate", "values"]
    for (k, row) in ((:liqprob, "liqprob"), (:d2c, "b2a"), (:CFrel, "CF reliance"))
        SH[Symbol(k, :_SME)][ind] = eq.sme[row, "SME Values"]
        SH[Symbol(k, :_LE)][ind]  = eq.sme[row, "LE Values"]
    end
end
shbar(k, title) = bar(shift_vec, SH[k], title = title, xlabel = "Shift Vec", legend = false)
savemisc(plot(shbar(:prodimp, "ProdImp"), shbar(:liqprob, "LiqProb"), shbar(:d2c, "D2C"), shbar(:CFrel, "CFRel"), shbar(:intrate, "IntRate"),
              shbar(:liqprob_SME, "LiqProb_SME"), shbar(:liqprob_LE, "LiqProb_LE"), shbar(:d2c_SME, "D2C_SME"), shbar(:d2c_LE, "D2C_LE"),
              shbar(:CFrel_SME, "CFRel_SME"), shbar(:CFrel_LE, "CFRel_LE"), layout = (4, 3), size = (1200, 900)), "M11_theta_shift")


########## Not ported ##########
# - ShiftingFixedCosts.jl: the lender pays a fixed cost in reorganization (Pi_reo = max(phi_c*V - zeta_Rs, 0)); a model change.
# - DebtScedule (M2): follows an older model version.
# - FindZeta.jl and SmallCosts.jl's grid search for zeta_Rs: replaced by find_zeta_Rs (Functions/Results.jl).
# - dynsim.jl: included by baseline_calib.jl but never called.
# - The calibration with ErrorFunc1 (Baseline.jl lines 655-684): replaced by calibrate_lazy (Functions/FCmodel/Calibrate.jl).

reset_params!(); set_params!(BaseEq.values)
println("MiscResults figures saved to ", misc_dir)

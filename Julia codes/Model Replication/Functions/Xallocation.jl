# FIRM AND CAPITAL ALLOCATION ACROSS PRODUCTIVITY STATES (fig:firmalloc and fig:capalloc in main5)
# The results block of FinalVers/Xallocation.jl as a function. That script solved its 4 scenarios itself, at fixed wages;
# here they come solved, from solve_eq or ge_eq (Functions/Results.jl):
#   runs: label => equilibrium for the baseline, reform (ζ/2) and reform (ζ/10), in this order
#   perf: label => the perfect credit equilibrium (perf_cred = 1), shown in the firm distribution panels only.
#         The capital panel shows the efficient allocation k_share_opt as "perfect credit", as the script did.
# Returns the capital allocation (plt1) and firmalloc, the mass (plt3) and share (plt2) of firms by productivity side
# by side as in main5's figure.
function Xallocation(runs, perf)

    # model sizes and parameters
    B = last(runs[1])
    x_size, e_size, k_size, b_size = B.gridsize
    rho_e, sigma_e, nul_e, alpha, nu, pc, beta, delta, pdef_exo_s, pdef_exo_l, discount, phi_a, phi_af, tauchen_sd, kappa, zeta_L, tau_vec, phi_c_hh, Fcut = B.params

    # productivity grid and "optimal" allocation objects (common across scenarios)
    prod_toskip = 0
    toskip = x_size * prod_toskip + 1

    e_chain = tauchen(e_size, rho_e, sigma_e, (1 - rho_e) * nul_e, tauchen_sd)
    e_vals  = exp.(e_chain.state_values)
    x       = e_vals[prod_toskip+1:end]

    e_stat_full = sum(stationary_distributions(e_chain))
    e_stat      = e_stat_full[prod_toskip+1:end]
    e_stat     /= sum(e_stat)

    k_share_opt = e_vals[prod_toskip+1:end].^(1/(1 - alpha - nu))
    k_share_opt ./= sum(k_share_opt)
    k_share_opt = e_stat .* k_share_opt
    k_share_opt ./= sum(k_share_opt)

    # ------------ Collect results (scenarios 1–4) ------------
    results = Dict{Int, Dict{Symbol,Any}}()
    for (id, (label, eq)) in enumerate([runs; perf])
        results[id] = Dict(
            :label => label,
            :SumPol => eq.SumPol,
            :mu => eq.mu,
            :m => eq.m,
            :xpol => eq.xpol
        )
    end

    # ------------ Per-scenario derived objects ------------
    for (id, R) in results
        mu = R[:mu]
        Span = toskip:(length(mu) - 1)

        mu_stat = sum(reshape(mu[Span], x_size, :), dims=1)[:]               # firm mass by productivity

        k_num = results[id][:SumPol][Span, 3] .* mu[Span]
        k_share_obs = k_num ./ sum(k_num)
        k_share_obs = sum(reshape(k_share_obs, x_size, :), dims=1)[:]         # by productivity

        R[:mu_stat] = mu_stat
        R[:k_share_obs] = k_share_obs
    end

    # ------------ Plot 1: Capital allocation (observed vs optimal) ------------
    # --- plt1: Capital allocation (skip id=4). id=2 dotted ---
    R1 = results[1]
    R2 = results[2]
    R3 = results[3]

    plt1 = plot(x, k_share_opt, label="perfect credit", lw=3,  color=:black, size=(900,600), margin = 5mm)
    plot!(plt1, x, R1[:k_share_obs], label=R1[:label], lw=3, color=:blue3)
    plot!(plt1, x, R2[:k_share_obs], label=R2[:label], lw=3, color=:red, linestyle=:dashdot, dashes=(1.0, 4.0))
    plot!(plt1, x, R3[:k_share_obs], label=R3[:label], linestyle=:dash, lw=3, color=:green)  # dark violet

    xlabel!(plt1, "Productivity")
    ylabel!(plt1, "Capital Share")
    #title!(plt1, "Firm Mass Across Productivity States")
        plot!(plt1, xtickfont=font(10), ytickfont=font(10), guidefont=font(12), legendfont=font(10))

    # ------------ Plot 2: Firm distribution by productivity ------------
    # --- plt2: Firm distribution (include id=4) ---
    R4 = results[4]

    plt2 = plot()

    plot!(plt2, x, R1[:mu_stat] / sum(R1[:mu_stat]), label=R1[:label], lw=2.5, color=:blue3)
    plot!(plt2, x, R2[:mu_stat] / sum(R2[:mu_stat]), label=R2[:label], lw=2.5, color=:red, linestyle=:dashdot, dashes=(1.0, 4.0))
    plot!(plt2, x, R3[:mu_stat] / sum(R3[:mu_stat]), label=R3[:label], lw=2.5, color=:green, linestyle=:dash)
    plot!(plt2, x, R4[:mu_stat] / sum(R4[:mu_stat]), label=R4[:label], lw=2.5, color=:black)  # violet for scenario 4

    xlabel!(plt2, "Productivity")
    ylabel!(plt2, "Share of Firm Mass")
    #title!(plt2, "Firm Mass Across Productivity States")
    plot!(plt2, xtickfont=font(10), ytickfont=font(10),
          guidefont=font(12), legendfont=font(10))

    # ------------ Plot 3: Firm distribution by productivity ------------
    plt3 = plot()

    plot!(plt3, x, R1[:mu_stat] , label=R1[:label], lw=2.5, color=:blue3)
    plot!(plt3, x, R2[:mu_stat], label=R2[:label], lw=2.5, color=:red, linestyle=:dashdot, dashes=(1.0, 4.0))
    plot!(plt3, x, R3[:mu_stat], label=R3[:label], lw=2.5, color=:green, linestyle=:dash)
    plot!(plt3, x, R4[:mu_stat], label=R4[:label], lw=2.5, color=:black)  # violet for scenario 4

    xlabel!(plt3, "Productivity")
    ylabel!(plt3, "Firm Mass")
    #title!(plt3, "Firm Mass Across Productivity States")
    plot!(plt3, xtickfont=font(10), ytickfont=font(10),
          guidefont=font(12), legendfont=font(10))

    firmalloc = plot(plt3, plt2, layout=(1,2), size=(1200, 500), margin=10mm)

    return (capital = plt1, firmalloc = firmalloc)
end

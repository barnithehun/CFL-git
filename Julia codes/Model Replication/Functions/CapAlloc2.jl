# RESULTS FOR MISALLOCATION DECOMPOSITION AND CAPITAL ALLOCATION
# SEE MORE DETAILED RESULTS IN CapAlloc2_old.jl
# grid and parameter definitions
x_size, e_size, k_size, b_size = gridsize()
rho_e, sigma_e, nul_e, alpha, nu, pc, beta, delta, pdef_exo_s, pdef_exo_l, discount, phi_a, phi_af, tauchen_sd, kappa,  zeta_L, tau_vec, phi_c_hh, Fcut = parameters()

prod_toskip = 0
toskip = x_size * prod_toskip + 1
Span = toskip:(length(mu) - 1)

## sdfsd
# EFFICIENT
e_chain = tauchen(e_size, rho_e, sigma_e, (1-rho_e)*nul_e, tauchen_sd)
e_vals = exp.(e_chain.state_values) 
e_active = e_vals[prod_toskip+1:end]

e_stat  = sum(stationary_distributions(e_chain))
e_stat = e_stat[prod_toskip+1:end] * (1/sum(e_stat[prod_toskip+1:end])) 

# Capital shares across productivities
k_share_opt = e_active .^ (1/(1-alpha-nu)) ./ sum(e_active .^ (1/(1-alpha-nu)))
k_share_opt = e_stat .* k_share_opt * (1/sum(e_stat .* k_share_opt))


# BASELINE
#phi_c =  0.45771403436007907
#zeta_Rl =  2795.820471809075
#@elapsed SumPol, e_chain, Fmat = FirmOptim(wage, phi_c, zeta_Rl)
#c_e, f0 = EntryValue(SumPol, e_chain) ;
#mu, m, xpol = stat_dist(SumPol, Fmat, f0);

#zeta_Rl = 1400
#@elapsed SumPol0, e_chain0, Fmat0 = FirmOptim(1.013, phi_c, zeta_Rl)
#c_e0, f00 = EntryValue(SumPol0, e_chain0); 
#mu0, m0, xpol0 = stat_dist(SumPol0, Fmat0, f00);


# sums 
mu_stat = sum(reshape(mu[Span], x_size, :), dims=1)[:]
k_share_obs = SumPol[Span,3] .* mu[Span] ./ sum(SumPol[Span,3] .* mu[Span])
k_share_obs = sum(reshape(k_share_obs, x_size, :), dims=1)[:]

mu_stat0 = sum(reshape(mu0[Span], x_size, :), dims=1)[:]
k_share_obs0 = SumPol0[Span,3] .* mu0[Span] ./ sum(SumPol0[Span,3] .* mu0[Span])
k_share_obs0 = sum(reshape(k_share_obs0, x_size, :), dims=1)[:]

# Plotting Capital Allocation
# This is the compbination of two effects: 
 #    - change in policies, firms that are affected the most hold more capital
 #    - leftward shift in the productivity distribution, more mass on the low end of the productivity distribution
x = e_vals[prod_toskip+1:end]
plot(x, k_share_obs,  label="Baseline", linewidth=2.5, color="#1f77b4")
plot!(x, k_share_obs0, label="Reform",  linewidth=2.5, color="#2ca02c")
plot!(x, k_share_opt,  label="Optimal", linewidth=2.5, color="#d62728", linestyle=:dash)
xlabel!("Productivity")
ylabel!("Capital Share")
plot!(xtickfont=font(10), ytickfont=font(10),
      guidefont=font(12), legendfont=font(10))

# savefig("C:\\Users\\szjud\\OneDrive\\Asztali gép\\EBCs\\CFL-git\\Latex codes\\Plots\\capital_allocation.png")

#---------------------------------------------------
#---------------------------------------------------
# MISALLOCATION Decomposition
# productivity calculated as a residual (baseline)
totK =  transpose(mu[1:end-1])*SumPol[1:end-1,3]
totL =  transpose(mu[1:end-1])*SumPol[1:end-1,10] # Ns = Nd
totY =  transpose(mu[1:end-1])*SumPol[1:end-1,11]
YtoL = totY/totL

A_obs = totY / (totK^alpha * totL^nu)
# mass of firms (baseline)
M = sum(mu[Span])
# composition effect (baseline)
C = sum((mu_stat ./ M) .* e_active .^ (1/(1-alpha-nu)))
# efficient policies given productivity and firm mass distribution (baseline)
A_star = M^(1-alpha-nu) * C^(1-alpha-nu)
# inefficiency factor (baseline)
Omega = A_obs / A_star
# reconstruct Y/L (baseline check)
(totK/totY)^(alpha/(1-alpha)) * (M^(1-alpha-nu)*C^(1-alpha-nu)*Omega)^(1/(1-alpha)) * totL^((alpha+nu-1)/(1-alpha))

# ------------------------------------------------------------------
# Reform (suffix 0)
totK0 =  transpose(mu0[1:end-1])*SumPol0[1:end-1,3]
totL0 =  transpose(mu0[1:end-1])*SumPol0[1:end-1,10] # Ns = Nd
totY0 =  transpose(mu0[1:end-1])*SumPol0[1:end-1,11]
YtoL0 = totY0/totL0

# productivity calculated as a residual (reform)
A_obs0 = totY0 / (totK0^alpha * totL0^nu)
# mass of firms (reform)
M0 = sum(mu0[Span])
# composition effect (reform)
C0 = sum((mu_stat0 ./ M0) .* e_active .^ (1/(1-alpha-nu)))
# efficient policies given productivity and firm mass distribution (reform)
A_star0 = M0^(1-alpha-nu) * C0^(1-alpha-nu)
# inefficiency factor (reform)
Omega0 = A_obs0 / A_star0
# reconstruct Y/L (reform check)
(totK0/totY0)^(alpha/(1-alpha)) * (M0^(1-alpha-nu)*C0^(1-alpha-nu)*Omega0)^(1/(1-alpha)) * totL0^((alpha+nu-1)/(1-alpha))


# Decomposition: 
total_effect = log(YtoL0) - log(YtoL)   # or just sum(effects)
CapDep = (alpha/(1-alpha)) * (log(totK0/totY0) - log(totK/totY)) 
ExtMarg = ((1-alpha-nu)/(1-alpha)) * (log(M0) - log(M)) - 0.003
IntMarg = ((1-alpha-nu)/(1-alpha)) * (log(C0) - log(C))
PolMarg = (1/(1-alpha))*(log(Omega0) - log(Omega)) + 0.003

effects = [CapDep, ExtMarg, IntMarg, PolMarg, total_effect]
labels  = ["Capital deepening", "Extensive margin", "Intensive margin", "Policy margin",  "Δ log(Y/L)"]

# make colors, last one red
colors = vcat(fill(:steelblue, length(effects)-1), :red)

bar(labels,   effects,  legend = false,
    ylabel = "Contribution to log(Y/L)",
    title = "Decomposition of Δ log(Y/L)",
    color = colors)



# You SHOULD GET CALCULATE THE PERFECT CREDIT ALLOCATION FROM A ZERO FRICTION MODEL VERSION, BECAUSE THAT TAKES INTO ACCOUNT THE ENTRY MARGIN
# Distribution of firms in productivity states 
plot(x, mu_stat / M , label = "baseline", lw = 2)  # THIS DOES NOT ALWAYS ADD UP TO 1, LOOK INTO IT
plot!(x, mu_stat0 / M0, label = "reform", lw = 2)
xlabel!("e values")
ylabel!("Value")
title!("e_stat and mu_stat across productivity states")


sum(mu_stat ./ M)
sum(mu_stat0 ./ M0)
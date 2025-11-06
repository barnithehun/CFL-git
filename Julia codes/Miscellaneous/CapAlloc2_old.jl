# Theoretical Best Allocation
# For some reasone the TFP measeure does not work with estat but the aggregates are good

# grid and parameter definitions
x_size, e_size, k_size, b_size = gridsize()
rho_e, sigma_e, nul_e, alpha, nu, pc, beta, delta, pdef_exo_s, pdef_exo_l, discount, phi_a, phi_af, tauchen_sd, kappa,  zeta_L, tau_vec, phi_c_hh, Fcut = parameters()

prod_toskip = 12
toskip = x_size * prod_toskip + 1
span = toskip:(length(mu) - 1)

mu_stat = sum(reshape(mu[span], x_size, :), dims=1)[:]
mu_stat0 = sum(reshape(mu0[span], x_size, :), dims=1)[:]


# Optimal capital shares across productivities
e_chain = tauchen(e_size, rho_e, sigma_e, (1-rho_e)*nul_e, tauchen_sd)
e_vals = exp.(e_chain.state_values) 
e_active = e_vals[prod_toskip+1:end]

e_stat  = sum(stationary_distributions(e_chain))
e_stat = e_stat[prod_toskip+1:end] * (1/sum(e_stat[prod_toskip+1:end])) 
k_share_opt = e_active .^ (1/(1-alpha-nu)) ./ sum(e_active .^ (1/(1-alpha-nu)))
k_share_opt = e_stat .* k_share_opt * (1/sum(e_stat .* k_share_opt))

mu_stat0 = sum(reshape(mu0[span], x_size, :), dims=1)[:]
TFP_eff = sum(e_stat .*  e_active .^ (1/(1-alpha-nu)))^(1-alpha-nu)

# Observed capital shares across productivities

#----------------------------------------------------------------------------------------------
# BASELINE

# TFP calc from firm policies
k_share_obs = SumPol[span,3] .* mu[span] ./ sum(SumPol[span,3] .* mu[span])
k_share_obs = sum(reshape(k_share_obs, x_size, :), dims=1)[:]

l_share_obs = SumPol[span,10] .* mu[span] ./ sum(SumPol[span,10] .* mu[span]) 
l_share_obs = sum(reshape(l_share_obs, x_size, :), dims=1)[:]
mu_stat = sum(reshape(mu[span], x_size, :), dims=1)[:]

TFP_obs = sum( e_active .* mu_stat .^ (1-alpha-nu) .* k_share_obs .^ alpha .*  l_share_obs .^ nu) 

# TFP calc from aggregates
totK =  transpose(mu[span])*SumPol[span,3]
totL =  transpose(mu[span])*SumPol[span,10] # Ns = Nd
totY =  transpose(mu[span])*SumPol[span,11]
YtoL = totY / totL

TFP_agg = totY / (totK^alpha*totL^nu)

A_obs = totY / (totK^alpha*totL^nu)

L_within =  A_obs / TFP_obs

A_obs^(1/(1-alpha)) * (totK/totY)^(alpha/(1-alpha)) * totL^(((alpha+nu-1)/(1-alpha)))


#----------------------------------------------------------------------------------------------
# REFORM
# TFP calc from firm policies

k_share_obs0 = SumPol0[span,3] .* mu0[span] ./ sum(SumPol0[span,3] .* mu0[span])
k_share_obs0 = sum(reshape(k_share_obs0, x_size, :), dims=1)[:]

l_share_obs0 = SumPol0[span,10] .* mu0[span] ./ sum(SumPol0[span,10] .* mu0[span]) 
l_share_obs0 = sum(reshape(l_share_obs0, x_size, :), dims=1)[:]

mu_stat0 = sum(reshape(mu0[span], x_size, :), dims=1)[:]
TFP_obs0 = sum( e_active .* mu_stat0 .^ (1-alpha-nu) .* k_share_obs0 .^ alpha .*  l_share_obs0 .^ nu) 


# TFP calc from aggregates
totK0 =  transpose(mu0[span])*SumPol0[span,3]
totL0 =  transpose(mu0[span])*SumPol0[span,10] # Ns = Nd
totY0 =  transpose(mu0[span])*SumPol0[span,11]
YtoL0 = totY0 / totL0

A_obs0 = totY0 / (totK0^alpha*totL0^nu)
L_within0 =  (totY0 / (totK0^alpha*totL0^nu)) / TFP_obs0


YtoL0/YtoL # difference in YtoL
A_obs0 / A_obs # difference in residual TFP
TFP_obs0 / TFP_obs # difference in disaggregated TFP
TFP0 = totY0 / (totK0^alpha*totL0^nu)


A_obs^(1/(1-alpha)) * (totK/totY)^(alpha/(1-alpha)) * totL^(((alpha+nu-1)/(1-alpha)))
A_obs0^(1/(1-alpha)) * (totK0/totY0)^(alpha/(1-alpha)) * totL0^(((alpha+nu-1)/(1-alpha)))

#--------------------------------------------------------------------------
# Decomposition in capital deepening and productivity
YtoL0 / YtoL
A_obs0^(1/(1-alpha)) / A_obs^(1/(1-alpha)) 
(totK0/totY0)^(alpha/(1-alpha)) / (totK/totY)^(alpha/(1-alpha)) 
totL^(((alpha+nu-1)/(1-alpha))) / totL0^(((alpha+nu-1)/(1-alpha)))

# Decomposition of A_obs
A_obs0 / A_obs
TFP_obs0 / TFP_obs 

sum( e_active .* (mu_stat0) .^ (1-alpha-nu) .* k_share_obs .^ alpha .*  l_share_obs .^ nu) / TFP_obs # Improvement do to different distribution of firms 
sum( e_active .* mu_stat .^ (1-alpha-nu) .* k_share_obs0 .^ alpha .*  l_share_obs0 .^ nu) / TFP_obs # Improvement due to difference in the policies accross bins 

# Improvement misallocation within productivity bins 
L_within0 / L_within

# The difference in observed productivity due to the above changes so... 
A_obs0 / A_obs # must be equal to 
(L_within0 / L_within) * (sum( e_active .* mu_stat0 .^ (1-alpha-nu) .* k_share_obs .^ alpha .*  l_share_obs .^ nu) / TFP_obs) * (sum( e_active .* mu_stat .^ (1-alpha-nu) .* k_share_obs0 .^ alpha .*  l_share_obs0 .^ nu) / TFP_obs)
# which holds! 

# effect of the change of firm mass
(sum(mu0) / sum(mu)) ^ (1-alpha-nu)

# --------------------------------------------------------
# In logs, decomposing Y
log(YtoL0) - log(YtoL) 
(alpha / (1-alpha)) * (log(totK0/totY0) - log(totK/totY))
(1/(1-alpha)) * (log(A_obs0) - log(A_obs))

# In logs decomposing A_obs
log(A_obs0) - log(A_obs)
(1-alpha-nu) * (log(sum(mu0)) -  log(sum(mu)))  

#---------------------------------------------------
# Plotting Capital Allocation
x = 1:length(k_share_obs)   # or whatever your index is
plot(x, k_share_obs, label="Observed 1")
plot!(x, k_share_obs0, label="Observed 2")
plot!(x, k_share_opt, label="Optimal")


#--------------------------------------------------
# Within and across brackets misallocation
# These two lines do not match: 
TFP_obs = sum( e_active .* mu_stat .^ (1-alpha-nu) .* k_share_obs .^ alpha .*  l_share_obs .^ nu) 
TFP_agg = totY / (totK^alpha*totL^nu) 
TFP_agg / TFP_obs


# group by e
μe  = sum(reshape(mu[span], x_size, :), dims=1)[:] # mass of firms per e
Ke  = sum(reshape(SumPol[span,3]  .* mu[span], x_size, :), dims=1)[:] # total capital per e
Le  = sum(reshape(SumPol[span,10] .* mu[span], x_size, :), dims=1)[:] # total labor per e

Ye_exact = sum(reshape(SumPol[span,11] .* mu[span], x_size, :), dims=1)[:] # total output per e

# theoretical upper bound per bin from collapsed shares
Ye_bound = e_active .* (μe .^(1 - alpha - nu)) .* (Ke .^ alpha) .* (Le .^ nu)

# within-bin proportionality factor (≤ 1; =1 iff proportional inputs)
Ξe = Ye_exact ./ Ye_bound

# Aggregate ratios
ratio_baseline = sum(Ye_exact) / sum(Ye_bound)   # ≈ TFP / TFP_obs

# The reason is because in my calculation of TFP_obs, I assume that each firm withing the same productivity bin has the same k_i/K which is not true. Due to credit frictions, firms with the same e_i may have different capital-labor ratios depending on their cash on hand. 



# Iteration again: 

# productivity calculated as a residual (baseline)
A_obs
# mass of firms (baseline)
M = sum(mu)
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

# productivity calculated as a residual (reform)
A_obs0
# mass of firms (reform)
M0 = sum(mu0)
# composition effect (reform)
C0 = sum((mu_stat0 ./ M0) .* e_active .^ (1/(1-alpha-nu)))
# efficient policies given productivity and firm mass distribution (reform)
A_star0 = M0^(1-alpha-nu) * C0^(1-alpha-nu)
# inefficiency factor (reform)
Omega0 = A_obs0 / A_star0
# reconstruct Y/L (reform check)
(totK0/totY0)^(alpha/(1-alpha)) * (M0^(1-alpha-nu)*C0^(1-alpha-nu)*Omega0)^(1/(1-alpha)) * totL0^((alpha+nu-1)/(1-alpha))


# Decomposition: 
log(YtoL0) - log(YtoL)
CapDep = (alpha/(1-alpha)) * (log(totK0/totY0) - log(totK/totY)) 
ExtMarg = ((1-alpha-nu)/(1-alpha)) * (log(M0) - log(M))
IntMarg = ((1-alpha-nu)/(1-alpha)) * (log(C0) - log(C))
PolMarg = (1/(1-alpha))*(log(Omega0) - log(Omega))


# Plotting
effects = [CapDep, ExtMarg, IntMarg, PolMarg]
labels  = ["Capital deepening", "Extensive margin", "Intensive margin", "Policy margin"]

bar(
    fill(1, length(effects)),          # put all stacks at x = 1
    effects,
    group = labels,                    # one stack per effect
    bar_position = :stack,
    legend = :outerright,
    xticks = ([1], ["Δ log(Y/L)"]),
    ylabel = "Contribution to log(Y/L)"
)

# optional: dashed line showing total (should equal sum of effects)
hline!([sum(effects)], linestyle = :dash, label = "Sum = $(round(sum(effects), digits=3))")



###########################################
#### CODE OUTTAKES TO TRY IN ISOLATION ####
###########################################
using LinearAlgebra, Statistics, LaTeXStrings, Plots, QuantEcon, Roots, NamedArrays
using SparseArrays, Dates, XLSX, DataFrames, Distributions, Random

beta = 0.96
pdef = 0.1
gam  = 0
Pi_liq = 4082
Pi_reo = 38823
next_b = 10857.1
tau_vec = 0:0.1:1

# vectorized for efficiency
   q_tau = zeros(length(tau_vec))
    
    if next_b == 0   
        fill!(q_tau, beta)
    else
        q_tau .= (beta ./ next_b) .* ((1 .- pdef) .* next_b .+
            pdef .* min.(next_b, gam .* ((1 .- tau_vec) .* Pi_liq .+ tau_vec .* kappa .* Pi_liq) .+ 
            (1 .- gam) .* ((1 .- tau_vec) .* Pi_liq .+ tau_vec .* Pi_reo)))
    end

    q, tau_index = findmax(q_tau)

    if isapprox(q, maximum(q_tau), atol=0.001) # == won't work here, there will always be a small numerical diff
        tau = gam*Pi_reo > (1-gam)*Pi_liq ? 1 : 0
    else
        tau = tau_vec[tau_index]
    end

    return q, tau




fn_Tau_Q(pdef, gam, Pi_liq, Pi_reo, next_b, tau_vec)
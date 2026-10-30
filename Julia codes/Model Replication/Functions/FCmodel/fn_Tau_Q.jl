# q and tau for one state-action pair: q is the highest debt price over the CFL reliances in tau_vec,
# tau the reliance that attains it. beta and kappa are passed in by the caller (before, they were read
# from parameters() on every call), and no temporary arrays are built, since this runs once per (e, a) in every q iteration.
function fn_Tau_Q(pdef, gam, Pi_liq, Pi_reo, next_b, tau_vec, beta, kappa)

    q = 0.0; tau_index = 0; q_min = 0.0
    for (i, tau_i) in enumerate(tau_vec)
        if next_b == 0
            q_i = beta
        else
            q_i = (beta / next_b) * ((1 - pdef) * next_b +
                pdef * gam * min(next_b, ((1 - tau_i) * Pi_liq + tau_i * kappa * Pi_liq)) +
                pdef * (1 - gam) * min(next_b, ((1 - tau_i) * Pi_liq + tau_i * Pi_reo)))
        end

        # highest q, the first one on ties (as findmax), and the lowest q (as minimum)
        if i == 1 || isless(q, q_i)
            q, tau_index = q_i, i
        end
        q_min = i == 1 ? q_i : min(q_min, q_i)
    end

    if isapprox(q, q_min, atol=0.001) # == won't work here, there will always be a small numerical diff
        tau = (1-gam)*Pi_reo > gam*Pi_liq ? 1 : 0
    else
        tau = tau_vec[tau_index]
    end

    return q, tau
end

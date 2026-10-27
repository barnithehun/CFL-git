function fn_Tau_Q(pdef, gam, Pi_liq, Pi_reo, next_b, tau_vec)

    beta = parameters().beta
    kappa = parameters().kappa

    # vectorized for efficiency
    q_tau = zeros(length(tau_vec))
    if next_b == 0   
        fill!(q_tau, beta)
    else
        q_tau .= (beta ./ next_b) .* ((1 .- pdef) .* next_b .+
            pdef .* gam .* min.(next_b, ((1 .- tau_vec) .* Pi_liq .+ tau_vec .* kappa .* Pi_liq)) .+ 
            pdef .* (1 .- gam) .* min.(next_b, ((1 .- tau_vec) .* Pi_liq .+ tau_vec .* Pi_reo)))
    end

    q, tau_index = findmax(q_tau)

    if isapprox(q, minimum(q_tau), atol=0.001) # == won't work here, there will always be a small numerical diff
        tau = (1-gam)*Pi_reo > gam*Pi_liq ? 1 : 0
    else
        tau = tau_vec[tau_index]
    end

    return q, tau
end
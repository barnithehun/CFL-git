function UpdateLoop(values::Vector{Float64}, policies::Vector{Int64}, phi_c)
    
    rho_e, sigma_e, nul_e, alpha, nu, pc, beta, delta, pdef_exo, discount, phi_a, tauchen_sd, kappa, zeta_R, zeta_L, tau_vec = parameters()
   
    x_size, e_size, _ , _ = gridsize()
    n,m, e_vals,e_ptrans, x_grid, _,s_i_vals, a_vals,a_i_vals = GridMake()
    fn_L(k,e) =  (nu*e*k^alpha / wage)^(1/(1-nu))
    fn_Y(k, e) =  e*k^alpha*fn_L(k,e)^nu
    fn_Pi(k, e) = (1-nu)*fn_Y(k,e)-pc
    fn_X(k,b,e) =  fn_Pi(k, e) + (1-delta) * k - b
    fn_D(next_k, next_b, x, q) =  x - next_k + q * next_b
    fn_Gam(k,val) = Int( (phi_a*(1-delta)*k - zeta_L) >= (phi_c*val - zeta_R))

    for s_i in 1:n

        e_i = s_i_vals[s_i, 2]

        for a_i in 1:m
            # policies given (x,e)
            next_k = a_vals[a_i,1]
            next_b = a_vals[a_i,2]

            pdef = 0
            gam = 0
            Pi_reo = 0
            Pi_liq = max(0, phi_a*(1-delta)*next_k - zeta_L)

            for next_e_i in 1:e_size

                p_trans = e_ptrans[e_i, next_e_i]
                x_next = fn_X(next_k,next_b,e_vals[next_e_i])

                x_close = argmin(abs.(x_next .- x_grid))   
                xe_close = x_close + (next_e_i-1)*x_size
                
                next_def_close = a_i_vals[policies[xe_close], 3]
                val_close = values[xe_close]

                if x_next < x_grid[end] && x_next > x_grid[1]
                    
                    x_far = x_next > x_grid[x_close] ? x_close + 1 : x_close - 1
                    # finding the correspoing indicies     
                    xe_far = x_far + (next_e_i-1)*x_size

                    next_def_far = a_i_vals[policies[xe_far], 3]
                    val_far = values[xe_far]
                    
                    close_weight = abs(x_next - x_grid[x_far]) / (abs(x_next - x_grid[x_close]) + abs(x_next - x_grid[x_far]))
                    
                    # value needed only for Pi_reo and gam
                    val = close_weight*val_close + (1-close_weight)*val_far    

                    pdef += p_trans*(close_weight*next_def_close + (1-close_weight)*next_def_far)
                    gam += p_trans * fn_Gam(next_k, val) 
                    Pi_reo += p_trans * max(phi_c*val - zeta_R, 0)

                else # close_weight = 1
                    val = val_close
                    pdef += p_trans * next_def_close                  
                    gam += p_trans * fn_Gam(next_k, val)
                    Pi_reo += p_trans * max(phi_c*val - zeta_R, 0)                        
                end  
            
            end 

            # saving results for summary
            q, tau = fn_Tau_Q(pdef, gam, Pi_liq, Pi_reo, next_b, tau_vec)
            pdef_sa[s_i, a_i] = pdef
            q_sa[s_i,a_i] = q
            tau_sa[s_i,a_i] = tau
            gam_sa[s_i, a_i] = gam
            Pi_liq_sa[s_i, a_i] = Pi_liq
            Pi_reo_sa[s_i, a_i] = Pi_reo
        end
    end

    return  pdef_sa, q_sa, tau_sa, gam_sa, Pi_liq_sa, Pi_reo_sa 

end
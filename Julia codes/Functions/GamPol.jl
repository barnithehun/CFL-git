#### hypotheticals ####
# interest rate and debt financing strategy across different liquidation probability
function GamPol(SumPol, x_i, e_i)

    gam_vec = 0:0.1:1
    tau_vec = 0:1.0:1
    x_size = gridsize().x_size    
    xe_i = x_i + (e_i-1)*x_size

    x_val = Int(round(SumPol[xe_i, 1],digits = 0))
    e_val = round(SumPol[xe_i, 2],digits = 1)
    b_val = Int(round(SumPol[xe_i, 4],digits = 0))

    next_b = SumPol[xe_i, 4]
    pdef = SumPol[xe_i, 8]
    Pi_liq = SumPol[xe_i, 15]
    Pi_reo = SumPol[xe_i, 16]

    qCgam = zeros(length(gam_vec))
    tauCgam = zeros(length(gam_vec))

    for (ind, gam) in enumerate(gam_vec)
        q, tau = fn_Tau_Q(pdef, gam, Pi_liq, Pi_reo, next_b, tau_vec)
        qCgam[ind] = q
        tauCgam[ind] = tau
    end

    plot(title= "TFP = $e_val; CoH = $x_val, Debt = $b_val", titlefont=font(12), legend=:topright, xlabel="Liquidation probability", ylabel="Inverse interest rate",  ylim=(parameters().beta-0.022, parameters().beta))
    # Plot qCgam with dynamic line style based on tauCgam value
    for i in 2:length(gam_vec)

        if tauCgam[i] == 1
            plot!(gam_vec[i-1:i], qCgam[i-1:i], label="", color=:blue, linestyle=:solid, lw=3)
        elseif tauCgam[i] == 0
            plot!(gam_vec[i-1:i], qCgam[i-1:i], label="", color=:red, linestyle=:solid, lw=3)
        end
    end

    plot!([], [], color=:blue, linestyle=:solid, label="CF borrowing", lw=3)
    plot!([], [], color=:red, linestyle=:solid, label="AB borrowing", lw=3)

end

# interest rate, debt financing strategy and probability of default across leverage
function DebtScedule(SumPol, x_i, e_i; phi_c)

    # starting functions 
    fn_L(k,e) =  (nu*e*k^alpha / wage)^(1/(1-nu))
    fn_Y(k, e) =  e*k^alpha*fn_L(k,e)^nu
    fn_Pi(k, e) = (1-nu)*fn_Y(k,e)-pc
    fn_X(k,b,e) =  fn_Pi(k, e) + (1-delta) * k - b
    fn_Gam(val, zeta_R) = Int( 0 >= ((1-phi_c_hh)*val - zeta_R)) # optimal liquidation

    # parameters and grid
    _,_, e_vals,_, e_ptrans, x_grid, _, s_i_vals, a_vals,a_i_vals = GridMake()
    rho_e, sigma_e, nul_e, alpha, nu, pc, beta, delta, pdef_exo_s, pdef_exo_l, discount, phi_a, tauchen_sd, kappa, zeta_Rs, zeta_Rl, zeta_L, tau_vec, phi_c_hh, Fcut = parameters()

    # calling grid size
    x_size, e_size, k_size, b_size = gridsize()


    # setting up the grid 
    x_size = gridsize().x_size    
    xe_i = x_i + (e_i-1)*x_size    
    
    # for later, to label the plot
    x_val = Int(round(SumPol[xe_i, 1],digits = 0))
    e_val = round(SumPol[xe_i, 2],digits = 3)
    
    # setting up k
    next_k = SumPol[xe_i, 3]
    k_val = Int(round(next_k, digits = 0))

    # taking the possible values of debt 
    b_vec  = 0:100:8*next_k
    n = length(b_vec)


    q_a = zeros(n)
    pdef_a = zeros(n)
    gam_a = zeros(n)
    Pi_liq_a = zeros(n)
    Pi_reo_a = zeros(n)
    tau_a = zeros(n)
    levrg = zeros(n)
    # Probability of default, liquidation, PIliq and PIreo and implied q given optimal k', b' in each state
    for (ind, next_b) in enumerate(b_vec)
        # policies given (x,e)

        zeta_R = next_k >= Fcut ? zeta_Rl : zeta_Rs
        pdef_exo = next_k >= Fcut ? pdef_exo_l : pdef_exo_s

        pdef = 0
        gam = 0
        Pi_reo = 0
        Pi_liq = max(0, phi_a*(1-delta)*next_k - zeta_L)

        for next_e_i in 1:e_size

            p_trans = e_ptrans[e_i, next_e_i]
            x_next = fn_X(next_k, next_b, e_vals[next_e_i])

            x_close = argmin(abs.(x_next .- x_grid))   
            xe_close = x_close + (next_e_i-1)*x_size
            
            next_def_close = SumPol[xe_close, 7]
            val_close = SumPol[xe_close, end]

            if x_next < x_grid[end] && x_next > x_grid[1]
                
                x_far = x_next > x_grid[x_close] ? x_close + 1 : x_close - 1
                # finding the correspoing indicies     
                xe_far = x_far + (next_e_i-1)*x_size

                next_def_far = SumPol[xe_far, 7]
                val_far = SumPol[xe_close, end]
                
                close_weight = abs(x_next - x_grid[x_far]) / (abs(x_next - x_grid[x_close]) + abs(x_next - x_grid[x_far]))
                
                # value needed only for Pi_reo and gam
                val = close_weight*val_close + (1-close_weight)*val_far    

                pdef += p_trans*(close_weight*next_def_close + (1-close_weight)*next_def_far)
                gam += p_trans *  fn_Gam(val, zeta_R) 
                Pi_reo += p_trans * max(phi_c*val - zeta_R, 0)

            else # close_weight = 1
                val = val_close
                pdef += p_trans * next_def_close                  
                gam += p_trans *   fn_Gam(val, zeta_R) 
                Pi_reo += p_trans * max(phi_c*val - zeta_R, 0)                        
            end  
        
        end 
        q, tau = fn_Tau_Q(pdef, gam, Pi_liq, Pi_reo, next_b, tau_vec)
        
        # saving results for summary
        pdef_a[ind] = pdef
        q_a[ind] =  q
        tau_a[ind] = tau
        gam_a[ind] = gam
        Pi_liq_a[ind] = Pi_liq
        Pi_reo_a[ind] = Pi_reo
        levrg[ind] =  next_b / next_k
    end

    p1 = plot(title= "TFP = $e_val", titlefont=font(12), legend=:topright, xlabel="Leverage", ylabel="Values")
    plot!(levrg, q_a, label="q", color=:blue, lw=3)
    plot!(levrg, pdef_a, label="pdef", color=:red, linestyle=:dash, lw=3)
    
    p2 = plot(title= "TFP = $e_val", titlefont=font(12), legend=:topright, xlabel="Leverage", ylabel="Values")
    plot!(levrg, gam_a, label="gam", color=:blue, lw=3)
    plot!(levrg, tau_a, label="tau", color=:red, linestyle=:dash, lw=3)
    

    return p1, p2
    
end

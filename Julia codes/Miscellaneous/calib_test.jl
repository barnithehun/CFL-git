###########################################################################
####### CALIBRATION - UPDATE FIRM OPTIM TO YOUR LATEST VERSION !!! ########
###########################################################################

# e_vec
pvec = [7 11 13 19 23 27 33] 

pvec = 0.1:0.3:1
resultsTOT = zeros(8,length(pvec));
resultsSME = zeros(7,length(pvec));
resultsLE = zeros(6,length(pvec));
@elapsed for (iter_a, PVAL) in enumerate(pvec)

 println("Iteration: $iter_a, Value: $PVAL")

 function gridsize()
    # grid sizes - x, k, b should be even numbers!!
    x_size::Int = 46
    e_size::Int = 23
    k_size::Int = 32
    b_size::Int = 32
    return (x_size = x_size, e_size = e_size, k_size = k_size, b_size = b_size)
end

    wage = 1;
    phi_c = 0.15;
    zeta_Rl = 2000;

function parameters()
    rho_e::Float64 = 0.969
    sigma_e::Float64 = 0.146
    nul_e::Int = 1
    DRS::Float64 = 0.75
    alpha::Float64 = 1/3 * DRS
    nu::Float64 = 2/3 * DRS
    pc::Float64 =   20
    beta::Float64 = 0.96
    delta::Float64 = 0.05
    pdef_exo_l::Float64 =  0.015
    pdef_exo_s::Float64 =  pdef_exo_l
    discount::Float64 = beta
    phi_a::Float64 =  0.4
    tauchen_sd::Float64 = 3

    kappa::Float64 = PVAL       # capital recovery rate of CFL debt
    # zeta_Rl::Float64 =   4000
    zeta_Rs::Float64 =  zeta_Rl
    tau_vec::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}} = 0:1.0:1  # vector of CFL reliances
    zeta_L::Float64 = 0
    phi_c_hh::Float64 =  0.4073866012405985
    Fcut::Float64 = 250

    return (rho_e = rho_e, sigma_e = sigma_e, nul_e = nul_e, alpha = alpha,
            nu = nu, pc = pc, beta = beta, delta = delta, pdef_exo_s = pdef_exo_s, pdef_exo_l = pdef_exo_l,
            discount = discount, phi_a = phi_a, tauchen_sd = tauchen_sd,
            kappa = kappa, zeta_Rs = zeta_Rs, zeta_L = zeta_L, tau_vec = tau_vec,
             phi_c_hh = phi_c_hh, Fcut = Fcut)
end
####### FIRM OPTIM #######
function FirmOptim(wage, phi_c, zeta_Rl)

    rho_e, sigma_e, nul_e, alpha, nu, pc, beta, delta, pdef_exo_s, pdef_exo_l, discount, phi_a, tauchen_sd, kappa, zeta_Rs, zeta_L, tau_vec, phi_c_hh, Fcut = parameters()

    # calling grid size
    x_size, e_size, k_size, b_size = gridsize()

    # setting optimization functions
    fn_L(k,e) =  (nu*e*k^alpha / wage)^(1/(1-nu))
    fn_Y(k, e) =  e*k^alpha*fn_L(k,e)^nu
    fn_Pi(k, e) = (1-nu)*fn_Y(k,e)-pc
    fn_X(k,b,e) =  fn_Pi(k, e) + (1-delta) * k - b
    fn_D(next_k, next_b, x, q) =  x - next_k + q * next_b
    
    # Liquidation policies
    fn_Gam(val, zeta_R) = Int( 0 >= ((1-phi_c)*val - zeta_R)) # optimal liquidation

    # calling grid size
    x_size, e_size, k_size, b_size = gridsize()
    
    # Setting the state-space
    # productivity process
    e_chain = tauchen(e_size, rho_e, sigma_e, (1-rho_e)*nul_e, tauchen_sd)
    e_vals = exp.(e_chain.state_values) 
    e_ptrans = e_chain.p 

    # Log-grids
    k_grid = [0;exp.(range(log(10), log(10^5), k_size-1))]
    b_grid = [0;exp.(range(log(10), log(10^5), b_size-1))]  # no savings

    x_low = fn_X(k_grid[1],b_grid[end], e_vals[1])
    x_high = fn_X(k_grid[end], b_grid[1], e_vals[end])
    x_grid_low = sort(-exp.(range(log(10), log(-x_low), ceil(div(x_size, 3))))[2:end], rev = false) # set the negative part of the x-grid size to 1/3rd
    x_grid_high = exp.(range(log(10), log(x_high), (x_size - length(x_grid_low)-1))) 
    x_grid = [ x_grid_low; 0; x_grid_high ]

    # Define Q(n,m,n) matrix
    n = x_size * e_size  # all possible states
    m = k_size * b_size  # all possible actions

    # total number of possible states, +1 state for being quit 
    s_vals = [gridmake(x_grid, e_vals) zeros(n)]           
    s_vals = [s_vals; [0 0 1]]

    s_i_vals = [gridmake(1:x_size, 1:e_size) zeros(n)] 
    s_i_vals = Int.([s_i_vals; [div(x_size,2) 1 1]])  # productivity after default does not matter

    a_vals = [gridmake(k_grid, b_grid) zeros(m,2) ]
    a_vals = [a_vals; [0 0 1 0]; [0 0 0 1]]

    a_i_vals = [gridmake(1:k_size, 1:b_size) zeros(m,2)]
    a_i_vals = Int.([a_i_vals; [1 div(b_size,2) 1 0]; [1 div(b_size,2) 0 1]])

    # adjusting the gridsize
    n = n+1
    m = m+2

    #################################### 
    # Collecting the probabilities of endogenous default in a matrix
    Pdefmat = zeros(n, m);
    for a_i in 1:m
        for  s_i in 1:n 
            next_k = a_vals[a_i, 1]   
            Pdefmat[s_i, a_i] = next_k >= Fcut ? pdef_exo_l : pdef_exo_s
        end
    end


    Q = zeros(n, m, n); 
    # exogenous default - no matter the state or action you will have a P_exo chance to end up in default the next period
    Q[:,:,end] = Pdefmat
    for a_i in 1:m
        for  s_i in 1:n 

            # productivities (indices)
            e = s_i_vals[s_i, 2]  # enough to save e, current x does not matter, since there are no financial frictions

            # actions (values)
            next_def = a_vals[a_i, 3] 
            next_exit = a_vals[a_i, 4] 
            def = s_vals[s_i,3]
            b = a_vals[a_i, 2]   
            k = a_vals[a_i, 1]   

            for next_e_i in 1:e_size  
                if  def == 0 && next_def == 0 && next_exit == 0
                
                    # next period cash on hand x'(b',k',e')
                    x_next = fn_X(k,b,e_vals[next_e_i])
                    # where x falls on the grid - closer
                    x_close = argmin(abs.(x_next .- x_grid))
            
                    # probability of transition from e_i to next_e_j - conditional on no def
                    pdef_exo = k >= Fcut ? pdef_exo_l : pdef_exo_s
                    p_trans = e_ptrans[e, next_e_i] * (1-pdef_exo)
            
                    # find the second closest
                    if x_next < x_grid[end] && x_next > x_grid[1]
                        
                        x_far = x_next > x_grid[x_close] ? x_close + 1 : x_close - 1            
                        close_weight = abs(x_next - x_grid[x_far]) / (abs(x_next - x_grid[x_close]) + abs(x_next - x_grid[x_far]))
            
                        # finding the correspoing indicies     
                        xe_close = x_close + (next_e_i-1)*x_size
                        xe_far = x_far + (next_e_i-1)*x_size
                        
                        # filling the transition matrix
                        Q[s_i, a_i, xe_close] = p_trans*close_weight
                        Q[s_i, a_i, xe_far] = p_trans*(1-close_weight)

                    else
                        xe_close = x_close + (next_e_i-1)*x_size
                        Q[s_i, a_i, xe_close] = p_trans
                    end  
                        
                else
                    # this states two things: 
                    # 1) if you are in an def state, you will be in an def state in the next period no matter the action (row of ones)
                    # 2) if you choose def or exit action, you will be in an def state in the next period no matter the state (columns of ones)
                    Q[s_i, a_i, end] = 1
                            
                end
            end
        end
    end

    # initial (!) endogenous default probability for each state
    kbexq_old::Array{Float64, 2} = zeros(n, 4)
    kbexq_new::Array{Float64, 2} = fill(1.0, n, 4)
    SumPol = zeros(n, 18)
    q_sa::Array{Float64, 2} = zeros(n, m)
    pdef_sa::Array{Float64, 2} = zeros(n, m)
    gam_sa::Array{Float64, 2} = zeros(n, m)
    Pi_liq_sa::Array{Float64, 2} = zeros(n, m)
    Pi_reo_sa::Array{Float64, 2} = zeros(n, m)
    tau_sa::Array{Float64, 2} = zeros(n, m)
    iter::Int = 0
    ################ 
    while !isequal(kbexq_old,kbexq_new)

        iter += 1
        if iter >= 25
            println("Error: Iteration number exceeded $iter")
            break
        end

        kbexq_old = kbexq_new
        R = fill(-Inf,  n, m);           
        for a_i in 1:m       
         # actions
         next_b = a_vals[a_i,2]
         next_k = a_vals[a_i,1]
         next_def = a_i_vals[a_i,3]
         next_exit =  a_i_vals[a_i,4]

            for s_i in 1:n
                
                def = s_vals[s_i, 3]
                x = s_vals[s_i, 1]

                if next_def == 0 && def == 0 && next_exit == 0
                        
                    # dividends
                    q = q_sa[s_i,a_i]
                    d = fn_D(next_k, next_b, x, q)

                    if d >= 0
                        R[s_i, a_i] = d 
                    end

                elseif next_def == 1 
                    R[s_i, a_i] = -5 
                elseif next_exit == 1
                    R[s_i, a_i] = x  
                elseif def == 1      
                    R[s_i, a_i] = 0
                end
            end
        end

        ddp = QuantEcon.DiscreteDP(R, Q, discount);
        results = QuantEcon.solve(ddp, PFI)

        Values = results.v;
        policies = results.sigma;  # optimal policy     

        ###################################################################
        # summarising results
        for s_i in 1:n

            # states
            x = s_vals[s_i, 1]
            e = s_vals[s_i, 2]

            # policies
            pol = policies[s_i]

            k = a_vals[pol,1]
            b = a_vals[pol,2]
            def = a_vals[pol,3]
            exit = a_vals[pol,4]

            # cash on hand if productivity stays the same
            next_x = fn_X(k, b, e)

            # implied firm policies
            l = fn_L(k,e)  # n is taken by gridsizefun
            y = fn_Y(k,e)
            Pi = fn_Pi(k,e)

            # if the firm defaults or exits these cannot be interpreted
            if def == 0 && exit == 0
                q = q_sa[s_i, pol]
                tau = tau_sa[s_i, pol]
                gam = gam_sa[s_i, pol]    
                pdef = pdef_sa[s_i, pol]
                Pi_liq = Pi_liq_sa[s_i, pol] 
                Pi_reo = Pi_reo_sa[s_i, pol] 
            else
                q = gam = Pi_liq = Pi_reo = tau = pdef = 0 
            end

            # if the firm defaults these cannot be interpreted
            d = def == 0 ? fn_D(k, b, x, q) : 0

            # value
            val = Values[s_i]

            # Summarise policies
            SumPol[s_i, :] .= [x, e, k, b, next_x, exit, def, pdef, q, l, y, Pi, d, gam, Pi_liq, Pi_reo, tau, val]    
        end
        
        ###############################################################################
        # Probability of default, liquidation, PIliq and PIreo and implied q given optimal k', b' in each state
        for s_i in 1:n

            x_i = s_i_vals[s_i, 1]
            e_i = s_i_vals[s_i, 2]
    
            for a_i in 1:m
                # policies given (x,e)
                next_k = a_vals[a_i,1]
                next_b = a_vals[a_i,2]

                # exo probability of default and reorganization costs change with size
                zeta_R = next_k >= Fcut ? zeta_Rl : zeta_Rs
                pdef_exo = next_k >= Fcut ? pdef_exo_l : pdef_exo_s
    
                pdef_endo = 0
                gam = 0
                Pi_reo = 0
                Pi_liq = max(0, phi_a*(1-delta)*next_k - zeta_L)

                for next_e_i in 1:e_size
    
                    p_trans = e_ptrans[e_i, next_e_i]
                    x_next = fn_X(next_k, next_b, e_vals[next_e_i])
    
                    x_close = argmin(abs.(x_next .- x_grid))   
                    xe_close = x_close + (next_e_i-1)*x_size
                    
                    next_def_close = a_i_vals[policies[xe_close], 3]
                    val_close = Values[xe_close]
    
                    if x_next < x_grid[end] && x_next > x_grid[1]
                        
                        x_far = x_next > x_grid[x_close] ? x_close + 1 : x_close - 1
                        # finding the corresponding indices     
                        xe_far = x_far + (next_e_i-1)*x_size
    
                        next_def_far = a_i_vals[policies[xe_far], 3]
                        val_far = Values[xe_far]
                        
                        close_weight = abs(x_next - x_grid[x_far]) / (abs(x_next - x_grid[x_close]) + abs(x_next - x_grid[x_far]))
                        
                        # value needed only for Pi_reo and gam
                        val = close_weight*val_close + (1-close_weight)*val_far    
    
                        pdef_endo += p_trans*(close_weight*next_def_close + (1-close_weight)*next_def_far)
                        gam += p_trans * fn_Gam(val, zeta_R) 
                        Pi_reo += p_trans * phi_c*val
    
                    else # close_weight = 1
                        val = val_close
                        pdef_endo += p_trans * next_def_close                  
                        gam += p_trans * fn_Gam(val, zeta_R)
                        Pi_reo += p_trans * phi_c*val         
                    end  
                
                end

                # now default is a decision and a shock at the same time
                # prob. of either pdef of pdef exo occurs- these are independent events
                pdef = pdef_endo + pdef_exo - pdef_endo * pdef_exo

                # adjusting the Pi_reo to the possibility of exo. default shock
                Pi_reo = Pi_reo * (1-pdef_exo)

                # updating Gamma: probability of liquidation either through endogenous or exogenous
                gam = (gam*pdef_exo + pdef_endo - ( gam*pdef_exo * pdef_endo)) / pdef

                # q and tau
                q, tau = fn_Tau_Q(pdef, gam, Pi_liq, Pi_reo, next_b, tau_vec)

                # saving results for summary
                pdef_sa[s_i, a_i] = pdef
                q_sa[s_i,a_i] = q
                tau_sa[s_i,a_i] = tau
                gam_sa[s_i, a_i] = gam
                Pi_liq_sa[s_i, a_i] = Pi_liq
                Pi_reo_sa[s_i, a_i] = Pi_reo
            end
        end
        kbexq_new = SumPol[:, [3, 4, 7, 9]]

    end
    println("Total 'main loop' iterations: ", iter)

    ### Incumbent dynamics ### 
    # Fmat - from state n, what is the probability of ending up in state n', given optimal policy
    Fmat = zeros(n-1,n-1)
    Fmat_bottom = zeros(n-1)
    for s_i in 1:(n-1)
               
        # policies imported from SumPol
        next_k = SumPol[s_i, 3]
        next_b = SumPol[s_i, 4]

        pdef_exo = next_k >= Fcut ? pdef_exo_l : pdef_exo_s

        Fmat_bottom[s_i] = pdef_exo + SumPol[s_i, 8]

        e_i = Int(floor( (s_i-1) / x_size) + 1) 
        for next_e_i in 1:e_size

            p_trans = e_ptrans[e_i, next_e_i] * (1-(pdef_exo+SumPol[s_i, 8]))
            x_next = fn_X(next_k,next_b,e_vals[next_e_i])

            x_close = argmin(abs.(x_next .- x_grid))
            xe_close = x_close + (next_e_i-1)*x_size

            if x_next < x_grid[end] && x_next > x_grid[1]
                
                x_far = x_next > x_grid[x_close] ? x_close + 1 : x_close - 1
                
                close_weight = abs(x_next - x_grid[x_far]) / (abs(x_next - x_grid[x_close]) + abs(x_next - x_grid[x_far]))
    
                # finding xe_index     
                xe_far = x_far + (next_e_i-1)*x_size

                # filling Fmat
                Fmat[s_i, xe_close] += p_trans * close_weight
                Fmat[s_i, xe_far] += p_trans * (1-close_weight)                               

            else # close_weight = 1
                Fmat[s_i, xe_close] += p_trans
            end  
        end

    end

    
    # Taking defaults into account - makes Fmat nXn
    Fmat = hcat(Fmat, Fmat_bottom)
    Fmat = vcat(Fmat,  [zeros(1,n-1) 1] )

   return ( SumPol, e_chain, transpose(Fmat) )

end
############ Results: calculation, only abl -> 0 ############ 
SumPol, e_chain, Fmat = FirmOptim(wage, phi_c, zeta_Rl)
c_e, f0 = EntryValue(SumPol, e_chain) 
mu, m, xpol = stat_dist(SumPol, Fmat, f0)

######## Saving results #######
n = size(SumPol,1)
totalmass = sum(mu)

# this contains firms that exit immidiately which are not counted in the model
exitmass=transpose(mu)*xpol
entrymass = m*(1-transpose(xpol)*f0) 
exitshare = exitmass/totalmass 

defshare = transpose(mu[1:n-1,:])*SumPol[1:n-1,7]/totalmass 

totK =  transpose(mu[1:n-1])*SumPol[1:n-1,3]
totB =  transpose(mu[1:n-1])*SumPol[1:n-1,4]
totL =  transpose(mu[1:n-1])*SumPol[1:n-1,10] # Ns = Nd
totY =  transpose(mu[1:n-1])*SumPol[1:n-1,11]

YtoL = totY/totL

# here LIE does not work bc. Im averaging ratios - loop is more readible than the vectorized version
avg_b2a, mu_b2a, avg_gam, mu_gam, avg_q, mu_q, avg_prod, mu_prod, avg_CFL, mu_CFL, share, CFdebt, Totdebt,  = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
for s_i in 1:n-1

    CFdebt += mu[s_i] * SumPol[s_i,4] * SumPol[s_i,17]
    Totdebt += mu[s_i] * SumPol[s_i,4]

    if SumPol[s_i,4] != 0
        avg_q += mu[s_i]/totalmass * SumPol[s_i,9]
        mu_q += mu[s_i]/totalmass
    end

    if SumPol[s_i,8] != 0 && SumPol[s_i,4] != 0
        avg_gam += mu[s_i]/totalmass * SumPol[s_i,14]
        mu_gam += mu[s_i]/totalmass
    end

    if SumPol[s_i,3] != 0 
        avg_b2a += mu[s_i]/totalmass * SumPol[s_i,4] / (SumPol[s_i, 3] + SumPol[s_i, 1])
        # avg_b2a += mu[s_i]/totalmass * SumPol[s_i,4] / (SumPol[s_i, 3] )
        mu_b2a += mu[s_i]/totalmass
    end  

    if SumPol[s_i,11] != 0 && SumPol[s_i, 10] != 0
        avg_prod += mu[s_i]/totalmass * SumPol[s_i,11] / SumPol[s_i, 10]
        mu_prod += mu[s_i]/totalmass
    end  

    if  SumPol[s_i, 4] != 0
        avg_CFL += mu[s_i]/totalmass * SumPol[s_i,17]
        mu_CFL += mu[s_i]/totalmass
    end  

    share += mu[s_i]/totalmass

    
end
avg_q = round((1/(avg_q / mu_q) - 1)*100, digits = 2)
avg_spread =  round(avg_q - (1/parameters().beta-1)*100, digits = 2)
avg_b2a = avg_b2a / mu_b2a
avg_gam = avg_gam / mu_gam
avg_prod = avg_prod / mu_prod # sanity check
avg_CFL = avg_CFL / mu_CFL
share = share / (1-exitshare) # mu_CFL has the same condition in the denominator
CFshare = CFdebt ./ Totdebt


resultsTOT[:,iter_a] = vcat(exitshare, defshare, avg_b2a, avg_q, avg_spread, avg_gam, CFshare, avg_CFL)

# here LIE does not work bc. Im averaging ratios - loop is more readible than the vectorized version
avg_b2a, mu_b2a, avg_gam, mu_gam, avg_q, mu_q, avg_prod, mu_prod, avg_CFL, mu_CFL, share, CFdebt, Totdebt,  = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
for s_i in 1:n-1

    if SumPol[s_i,3] < parameters().Fcut

        CFdebt += mu[s_i] * SumPol[s_i,4] * SumPol[s_i,17]
        Totdebt += mu[s_i] * SumPol[s_i,4]

        if SumPol[s_i,4] != 0
            avg_q += mu[s_i]/totalmass * SumPol[s_i,9]
            mu_q += mu[s_i]/totalmass
        end

        if SumPol[s_i,8] != 0 && SumPol[s_i,4] != 0
            avg_gam += mu[s_i]/totalmass * SumPol[s_i,14]
            mu_gam += mu[s_i]/totalmass
        end

        if SumPol[s_i,3] != 0 
            avg_b2a += mu[s_i]/totalmass * SumPol[s_i,4] / (SumPol[s_i, 3] + SumPol[s_i, 1])
            # avg_b2a += mu[s_i]/totalmass * SumPol[s_i,4] / (SumPol[s_i, 3] )
            mu_b2a += mu[s_i]/totalmass
        end  

        if SumPol[s_i,11] != 0 && SumPol[s_i, 10] != 0
            avg_prod += mu[s_i]/totalmass * SumPol[s_i,11] / SumPol[s_i, 10]
            mu_prod += mu[s_i]/totalmass
        end  

        if  SumPol[s_i, 4] != 0
            avg_CFL += mu[s_i]/totalmass * SumPol[s_i,17]
            mu_CFL += mu[s_i]/totalmass
        end  

        share += mu[s_i]/totalmass

    end
    
end
avg_q = round((1/(avg_q / mu_q) - 1)*100, digits = 2)
avg_spread =  round(avg_q - (1/parameters().beta-1)*100, digits = 2)
avg_b2a = avg_b2a / mu_b2a
avg_gam = avg_gam / mu_gam
avg_prod = avg_prod / mu_prod # sanity check
avg_CFL = avg_CFL / mu_CFL
share = share / (1-exitshare) # mu_CFL has the same condition in the denominator
CFshare = CFdebt ./ Totdebt

resultsSME[:,iter_a] = vcat(avg_b2a, avg_q, avg_spread, avg_gam, CFshare, avg_CFL, share )

# here LIE does not work bc. Im averaging ratios - loop is more readible than the vectorized version
avg_b2a, mu_b2a, avg_gam, mu_gam, avg_q, mu_q, avg_prod, mu_prod, avg_CFL, mu_CFL, share, CFdebt, Totdebt,  = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
for s_i in 1:n-1

    if SumPol[s_i,3] >= parameters().Fcut

        
        CFdebt += mu[s_i] * SumPol[s_i,4] * SumPol[s_i,17]
        Totdebt += mu[s_i] * SumPol[s_i,4]

        if SumPol[s_i,4] != 0
            avg_q += mu[s_i]/totalmass * SumPol[s_i,9]
            mu_q += mu[s_i]/totalmass
        end

        if SumPol[s_i,8] != 0 && SumPol[s_i,4] != 0
            avg_gam += mu[s_i]/totalmass * SumPol[s_i,14]
            mu_gam += mu[s_i]/totalmass
        end

        if SumPol[s_i,3] != 0 
            avg_b2a += mu[s_i]/totalmass * SumPol[s_i,4] / (SumPol[s_i, 3] + SumPol[s_i, 1])
            # avg_b2a += mu[s_i]/totalmass * SumPol[s_i,4] / (SumPol[s_i, 3] )
            mu_b2a += mu[s_i]/totalmass
        end  

        if SumPol[s_i,11] != 0 && SumPol[s_i, 10] != 0
            avg_prod += mu[s_i]/totalmass * SumPol[s_i,11] / SumPol[s_i, 10]
            mu_prod += mu[s_i]/totalmass
        end  

        if  SumPol[s_i, 4] != 0
            avg_CFL += mu[s_i]/totalmass * SumPol[s_i,17]
            mu_CFL += mu[s_i]/totalmass
        end  

        share += mu[s_i]/totalmass

    end
    
end
avg_q = round((1/(avg_q / mu_q) - 1)*100, digits = 2)
avg_spread =  round(avg_q - (1/parameters().beta-1)*100, digits = 2)
avg_b2a = avg_b2a / mu_b2a
avg_gam = avg_gam / mu_gam
avg_prod = avg_prod / mu_prod # sanity check
avg_CFL = avg_CFL / mu_CFL
share = share # mu_CFL has the same condition in the denominator
CFshare = CFdebt ./ Totdebt

resultsLE[:,iter_a] = vcat( avg_b2a, avg_q, avg_spread, avg_gam, CFshare, avg_CFL)

end

varnames = ["exitshare", "defshare", "TOT_b2a", "TOT_intrate", "TOT_spread",  "TOT_liqprob", "TOT_CFshare", "TOT_CFL", "LE_b2a", "LE_intrate", "LE_spread",  "LE_liqprob", "LE_CFshare", "LE_CFL", "SME_b2a", "SME_intrate", "SME_spread",  "SME_liqprob", "SME_CFshare", "SME_CFL", "SME_share"];

pvec_str =string.(vec(pvec));
results = vcat(resultsTOT, resultsLE , resultsSME);
results = NamedArray(results, names=( varnames,  pvec_str ) ,  dimnames=("Res", "ParamVal"))
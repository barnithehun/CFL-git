########################################################################
########### VER 7.4 - exogeneous liquidation probability ###############
########################################################################

using LinearAlgebra, Statistics, LaTeXStrings, Plots, QuantEcon, Roots, NamedArrays, Dates, XLSX, DataFrames, Distributions, Random, CSV

function gridsize()
    # grids sizes - x,k,b should be even numbers!!
    return(
        x_size = 50,
        e_size = 11,
        k_size = 40,
        b_size = 40)
end

function parameters()
    rho_e = 0.969
    sigma_e = 0.146
    nul_e = 1
    DRS = 0.75
    alpha = 1/3 * DRS
    nu = 2/3 * DRS
    pc = 1000
    beta = 0.96
    delta = 0.06
    pdef_exo = 0.04
    discount = beta
    phi_a = 0.4
    tauchen_sd = 4

    kappa = 0.3               # capital recovery rate of CFL debt
    zeta = 10^4               # fixed cost of reorganization
    tau_vec = 0:1             # vector of CFL reliances - optimal will be one or zero

    return (rho_e = rho_e, sigma_e = sigma_e, nul_e = nul_e, alpha = alpha,
            nu = nu, pc = pc, beta = beta, delta = delta, pdef_exo = pdef_exo,
            discount = discount, phi_a = phi_a, tauchen_sd = tauchen_sd, kappa = kappa, zeta = zeta, tau_vec = tau_vec)

end

####### FIRM OPTIM #######
function FirmOptim(wage; phi_c)

    # calling parameters
    rho_e, sigma_e, nul_e, alpha, nu, pc, beta, delta, pdef_exo, discount, phi_a, tauchen_sd, kappa, zeta, tau_vec = parameters()

    # calling grid size
    x_size, e_size, k_size, b_size = gridsize()

    # setting optimization functions
        fn_L(k,e) =  (nu*e*k^alpha / wage)^(1/(1-nu))
        fn_Y(k, e) =  e*k^alpha*fn_L(k,e)^nu
        fn_Pi(k, e) = (1-nu)*fn_Y(k,e)-pc
        fn_X(k,b,e) =  fn_Pi(k, e) + (1-delta) * k - b
        fn_D(next_k, next_b, x, q) =  x - next_k + q * next_b

        fn_chi(k,val) = Int(phi_a*(1-delta)*k >= phi_c*val - zeta) # liquidation decision b)

        function fn_Tau_Q(pdef, gam, Pi_liq, Pi_reo, next_b, tau_vec)

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
                tau = (1-gam)*Pi_reo > gam*Pi_liq ? 1 : 0
            else
                tau = tau_vec[tau_index]
            end

            return q, tau
        end

    # Setting the state-space
        # productivity process
        e_chain = tauchen(e_size, rho_e, sigma_e, (1-rho_e)*nul_e, tauchen_sd)
        e_vals = exp.(e_chain.state_values)
        # adding exogeneous default shocks
        e_ptrans = e_chain.p .* (1-pdef_exo) 
        e_ptrans[:,1] = e_ptrans[:,1] .+ pdef_exo

        # Log-grids
        k_grid = [0;exp.(range(log(10), log(10*10^5), k_size-1))]
        b_grid = [0;exp.(range(log(10), log(10*10^5), b_size-1))]  # no savings

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
        
    #################################### Standard approach
    Q = zeros(n, m, n); 
    for a_i in 1:m
        for  s_i in 1:n 

         # productivities (indicies)
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
            
                    # probability of transition from e_i to next_e_j 
                    p_trans = e_ptrans[e, next_e_i] 
            
                    # find the second closest
                    if x_next < x_grid[end] && x_next > x_grid[1]
                        
                        if x_next > x_grid[x_close]
                            x_far = x_close + 1
                            else
                            x_far = x_close - 1
                        end
                        
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
                    # 1) if you are in an def state, you will be in an def state in the next period no matter the action
                    # 2) if you choose and def action, you will be in an def state in the next period no matter the state
                    Q[s_i, a_i, end] = 1
                            
                end
            end
        end
    end

    # initital (!) endogeneous default probability for each state
    kbexq_old = zeros(n,4)
    kbexq_new = fill(1,n,4)
    SumPol = zeros(n, 18)

    q_sa = fill(0.97,n,m)
    pdef_sa = zeros(n,m);
    gam_sa = zeros(n,m);
    Pi_liq_sa = zeros(n,m);
    Pi_reo_sa =  zeros(n,m);
    tau_sa = zeros(n,m);
    iter = 0
    ################ 
    while !isequal(kbexq_old,kbexq_new)

        iter += 1
        if iter > 50
            println("Error: Iteration number exceeded $iter")
            break
        end

        kbexq_old = kbexq_new
        R = fill(-Inf,  n, m);           
        for a_i in 1:m       
         # actions
         next_b = a_vals[a_i,2]
         next_k = a_vals[a_i,1]
         next_def =  a_i_vals[a_i,3]
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
                    R[s_i, a_i] = -5  # -5000 if you want anyone to quit on its own
                elseif next_exit == 1
                    R[s_i, a_i] = x  
                elseif def == 1      
                    R[s_i, a_i] = 0
                end
            end
        end

        ddp = QuantEcon.DiscreteDP(R, Q, discount);
        results = QuantEcon.solve(ddp, PFI)

        values = results.v;
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
            pdef = pdef_sa[s_i, pol]

            # Def
            if def == 0
                q = q_sa[s_i, pol]
                tau = tau_sa[s_i, pol]
                gam = gam_sa[s_i, pol]    
                d = fn_D(k, b, x, q)
                Pi_liq = Pi_liq_sa[s_i, pol] 
                Pi_reo = Pi_reo_sa[s_i, pol] 
            else
                q = d = gam = Pi_liq = Pi_reo = tau = 0 
            end

            # value
            val = values[s_i]

            # Summarise policies
            SumPol[s_i, :] .= [x, e, k, b, next_x, exit, def, pdef, q, l, y, Pi, d, gam, Pi_liq, Pi_reo, tau, val]    
            
        end
        
        ###############################################################################
        # Probability of default, liquidation, PIliq and PIreo and implied q given optimal k', b' in each state
        for s_i in 1:n

            e_i = s_i_vals[s_i, 2]

            for a_i in 1:m
                # policies given (x,e)
                next_k = a_vals[a_i,1]
                next_b = a_vals[a_i,2]

                pdef = 0
                gam = 0
                Pi_reo = 0
                Pi_liq = phi_a*(1-delta)*next_k

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
                        gam += p_trans * fn_chi(next_k, val) 
                        Pi_reo += p_trans * max(phi_c*val - zeta, 0)

                    else # close_weight = 1
                        val = val_close
                        pdef += p_trans * next_def_close
                        gam += p_trans * fn_chi(next_k, val)
                        Pi_reo += p_trans * max(phi_c*val - zeta, 0)                        
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

        kbexq_new = SumPol[:, [3, 4, 7, 9]]

    end
    println("Total 'main loop' iterations: ", iter)

    ### Incumbent dynamics ### 
    # Fmat - from state n, what is the probability of ending up in state n', given optimal policy
    Fmat = zeros(n-1,n-1)
    for s_i in 1:(n-1)
               
        # policies imported from SumPol
        next_k = SumPol[s_i, 3]
        next_b = SumPol[s_i, 4]
        
        e_i = Int(floor( (s_i-1) / x_size) + 1) 
        for next_e_i in 1:e_size

            p_trans = e_ptrans[e_i, next_e_i]
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
    Fmat = Fmat .* (1 .- SumPol[1:end-1, 8])
    Fmat = hcat(Fmat, SumPol[1:end-1, 8])
    Fmat = vcat(Fmat,  [zeros(1,n-1) 1] )

   return ( SumPol, e_chain, transpose(Fmat) )

end

####### ENTRANTS #######
function EntryValue(SumPol, e_chain)  

    # entrant ln(e) distribution
    # here, set to be equal to the stationary distribution of e_chain    
    e_entry  = reduce(+,stationary_distributions(e_chain))

    # entrant X distribution - x_e = 0  in every case 
    x_vals = unique(SumPol[:, 1])
    zero_index = findall(x -> x == 0.0, x_vals)
    x_entry = zeros(length(x_vals))
    x_entry[zero_index .+ 0] .= 1 # if x = 0 prob = 1 

    # (x,e) are independent, the joint of the two distribution is their product
    xe_entry = [kron(e_entry, x_entry); 0] # also the f0 vector

    # map the entry probabilities to values
    beta = parameters().beta
    Ve = transpose(xe_entry) * (SumPol[:,end])*beta

    return ( Ve, xe_entry )    

end

wage = 1
@elapsed SumPol, e_chain, Fmat = FirmOptim(wage, phi_c = 0.7)
c_e, f0 = EntryValue(SumPol, e_chain) # free entry condition

####### STATIONARY DISTRIBUTION ########
function stat_dist(SumPol, Fmat, f0)

    # Exiting firms (voluntary + involuntary)
    n = size(SumPol,1)
    xpol = [SumPol[1:n-1,6] + SumPol[1:n-1,7] ; 1]
    Ident = Matrix(I,n,n)

    xpol_mat = Ident - Diagonal(xpol)   # I - diag(X)
    f0 = xpol_mat*f0                    # (I - diag(X))f0 - entrants may quit immidiately
    Mmat = Fmat*xpol_mat                # M = F(I - diag(X))

    # unscaled stationary distribution
    mu_0 = inv(Ident - Mmat)*f0         # inv(I-M)*f0

    # ok, bc exit and default implies k = 0, n = 0
    Nd = transpose(SumPol[1:n,10])*mu_0

    # This is given that labour supply is one (10000) inelastically
    m = 10000/Nd
    mu = m.*mu_0

    return ( mu, m , xpol )
        
end

# Print policies
function PrintPol()    
    mu, _, _ = stat_dist(SumPol, Fmat, f0)
    column_names = [:x, :e, :k, :b, :next_x, :exit, :def, :pdef, :q, :l, :y, :Pi, :d, :gam, :Pi_liq, :Pi_reo, :tau, :val, :SSpercent]
    SumPoll = hcat(SumPol[1:end, :],round.(mu ./ sum(mu) .* 100, digits=4))
    SumPol_df = DataFrame(SumPoll, column_names)
    time_string = "$(Dates.day(now()))_$(Dates.hour(now()))_$(Dates.minute(now()))_$(Dates.second(now()))"
    XLSX.writetable("$time_string.xlsx", sheetname="sheet1", SumPol_df)
end
# PrintPol()    

################ Results: steady state values ###################
    function sumSS(SumPol,Fmat,f0)
        
        mu, m, xpol = stat_dist(SumPol, Fmat, f0)
        n = size(SumPol,1)
        totalmass = sum(mu)

        # THIS IS PRLY NOT CORRECT BUT IT ADDS UP AT LEAST
        # this contains firms that exit immidiately which are not counted in the model
        exitmass=transpose(mu)*xpol
        entrymass = m*(1-transpose(xpol)*f0)  # m, needs to be adjusted with the firms that decide to quit immidiately
        exitshare = exitmass/totalmass 

        totK =  transpose(mu)*SumPol[1:n,3]
        totB =  transpose(mu)*SumPol[1:n,4]
        totL =  transpose(mu)*SumPol[1:n,10] # Ns = Nd
        totY =  transpose(mu)*SumPol[1:n,11]
        totPi = transpose(mu)*SumPol[1:n,12]

        YtoL = totY/totL

        meanL = totL / totalmass 
        meanK = totK / totalmass 
        meanY = totY / totalmass 

        # here LIE does not work bc. Im averaging ratios - loop is more readible than the vectorized version
        avg_b2k, mu_b2k, avg_q, mu_q, avg_prod, mu_prod, avg_CFL, mu_CFL, SMEshare = 0, 0, 0, 0, 0, 0, 0, 0, 0
        for s_i in 1:n

            if (SumPol[s_i,6] + SumPol[s_i,7]) != 1
                avg_q += mu[s_i]/totalmass * SumPol[s_i,9]
                mu_q += mu[s_i]/totalmass
            end

            if SumPol[s_i,4] != 0 && SumPol[s_i, 3] != 0
                avg_b2k += mu[s_i]/totalmass * SumPol[s_i,4] / SumPol[s_i, 3]
                mu_b2k += mu[s_i]/totalmass
            end  

            if SumPol[s_i,11] != 0 && SumPol[s_i, 10] != 0
                avg_prod += mu[s_i]/totalmass * SumPol[s_i,11] ./ SumPol[s_i, 10]
                mu_prod += mu[s_i]/totalmass
            end  

            if  SumPol[s_i, 3] != 0
                avg_CFL += mu[s_i]/totalmass * SumPol[s_i,17]
                mu_CFL += mu[s_i]/totalmass
            end  

            if SumPol[s_i, 3] != 0 && (SumPol[s_i, 3]+SumPol[s_i, 1]) <= 5000 # smake definition
                SMEshare += mu[s_i]/totalmass
            end
            
        end
        avg_q = avg_q / mu_q
        avg_b2k = avg_b2k / mu_b2k
        avg_prod = avg_prod / mu_prod # sanity check
        avg_CFL = avg_CFL / mu_CFL
        SMEshare = SMEshare / mu_CFL # mu_CFL has the same condition in the denominator

        CFshare = (transpose(mu)*(SumPol[1:n,4] .* SumPol[1:n,17]))  /  totB

        results = zeros(19,1)
        results[:,1] = vcat(totalmass, exitmass, entrymass, exitshare, totK, totB, totL, totY, totPi, YtoL, meanL, meanK, meanY, avg_b2k, avg_q, avg_prod, CFshare, avg_CFL, SMEshare)
        varnames = ["totalmass", "exitmass", "entrymass", "exitshare", "totK", "totB", "totL", "totY", "totPi", "YtoL", "meanL", "meanK", "meanY", "avg_b2k", "avg_q", "avg_prod", "CFshare", "avg_CFL", "SMEshare"];
        results = NamedArray(results, names=( varnames, ["values"] ) ,  dimnames=("Res", "ParamVal"))

    return ( results )

    end

    sumSS(SumPol,Fmat,f0)

############ Results: stationary distributions ############
    function StatDist(binnum, var::Char, SumPol)

        char_to_number = Dict('k' => 3, 'b' => 4, 'l' => 10, 'y' => 11, 'p' => 12, 'v' => 18)
        varnum = get(char_to_number, var, 0)
        max = maximum(SumPol[:, varnum])

        bins = [0; exp.(range(log(10), log(max+1), binnum-1))]
        binfill = zeros(binnum, 1)
        
        n = size(SumPol,1)-1
        for s_i = 1:n

            val = SumPol[s_i,varnum] 
            val_close = argmin(abs.(val .- bins))

            if val_close != binnum
                    
                if val > bins[val_close]
                    binfill[val_close + 1] += mu[s_i]
                else
                    binfill[val_close] += mu[s_i]
                end

            else
                binfill[val_close] += mu[s_i]
            end        
        end

        PDF = binfill ./ sum(binfill)
        CDF = cumsum(PDF, dims = 1)

        return ( PDF, CDF, bins ) 
        
    end
    mu, m , xpol = stat_dist(SumPol, Fmat, f0)


    function plotPDF(binnum, var::Char, SumPol)

        PDF, _, bins = StatDist(binnum, var::Char, SumPol)
        bin_labels = string.(round.(bins ./ 1000, digits=2))
        plot = bar(bin_labels, PDF, title=var, xrotation=45)
        return ( plot )

    end

    function plotCDF(binnum, var::Char, SumPol)

        _, CDF, bins = StatDist(binnum, var::Char, SumPol)
        bin_labels = string.(round.(bins ./ 1000, digits=2))
        plott = plot(bin_labels, CDF, title=var, xrotation=45, legend=false, linewidth=3)
        return ( plott )

    end

    binnum = 10
    plot(plotPDF(binnum, 'k', SumPol),
        plotPDF(binnum, 'b', SumPol),
        plotPDF(binnum, 'l', SumPol),
        plotPDF(binnum, 'y', SumPol),
        plotPDF(binnum, 'p', SumPol),
        plotPDF(binnum, 'v', SumPol), layout=(2,3), size=(1200, 800))
    plot(plotCDF(binnum, 'k', SumPol),
        plotCDF(binnum, 'b', SumPol),
        plotCDF(binnum, 'l', SumPol),
        plotCDF(binnum, 'y', SumPol),
        plotCDF(binnum, 'p', SumPol),
        plotCDF(binnum, 'v', SumPol), layout=(2,3), size=(1200, 800))


    # XE distribution
    function plotXE(SumPol, mu, e_chain)

        x_size, e_size, _, _ = gridsize()
            x_dist = zeros(x_size)
        for s_i in 1:e_size
            x_dist += mu[1 + (s_i-1)*x_size : s_i*x_size]
        end

        # stationary e distribution
        e_dist = zeros(e_size)
        for s_i in 1:e_size
            e_dist[s_i] = sum(mu[1 + (s_i-1)*x_size : s_i*x_size])
        end

        return (plot(bar(string.(round.(exp.(e_chain.state_values))), e_dist, title = "e_dist"),
                bar(string.(round.(unique(SumPol[:, 1])./1000)),  x_dist, title = "x_dist"), 
                layout=(2,1), size=(1200, 800)))

    end
    plotXE(SumPol, mu, e_chain)

############ Results: dynamics simulations ##############
    # In this version, there are productivity shocks, you could consider a setup absent of them
    function xkb_simul(SumPol, Fmat; simn_length = 10000, e_i)

        function next_si(Fvec)
            
            r = rand()  # uniform distribution between 0-1
            cumulative_prob = 0.0
            
            for (i, prob) in enumerate(Fvec)
                cumulative_prob += prob
                if r < cumulative_prob
                    return i
                end
            end 

            return length(Fvec)  # default to the last state if no transition
            
        end

        x_size, e_size, _, _ = gridsize()

        simt_length = 25
        simmat_k = fill(NaN, simn_length,simt_length);
        simmat_b = fill(NaN, simn_length,simt_length);
        simmat_x = fill(NaN, simn_length,simt_length);
        for simn in 1:simn_length 

            x_i = findall(x -> x == 0.0, unique(SumPol[:, 1]))[1]
            s_i = x_i + (e_i-1)*x_size

            for simt in 1:simt_length 

                xpol = SumPol[s_i,6] + SumPol[s_i,7] 

                if xpol != 1
                    
                    if simt == 1
                        simmat_x[simn,simt] = SumPol[s_i,1]
                        simmat_k[simn,simt] = SumPol[s_i,3]
                        simmat_b[simn,simt] = SumPol[s_i,4]
                    else
                        # Fmat is the transposed! transition matrix   
                        Fvec = Fmat[:,s_i]
                        s_i = next_si(Fvec)
                        simmat_x[simn,simt] = SumPol[s_i,1]
                        simmat_k[simn,simt] = SumPol[s_i,3]
                        simmat_b[simn,simt] = SumPol[s_i,4]
                    end

                else
                    break
                end
            end
        end

        meanX = zeros(simt_length)
        meanK = zeros(simt_length) 
        meanB = zeros(simt_length)
        n_share = zeros(simt_length)
        for simt in 1:simt_length 

            meanX[simt, 1] = mean(filter(!isnan, simmat_x[:, simt]))
            meanK[simt, 1] = mean(filter(!isnan, simmat_k[:, simt]))
            meanB[simt, 1] = mean(filter(!isnan, simmat_b[:, simt]))
            n_share[simt, 1] =   (simn_length - length(filter(!isnan, simmat_x[:, simt])))/simn_length

        end

        x_axis = 1:simt_length

        plott1 = plot(x_axis, meanX, label="X", linewidth=3)
        plot!(x_axis, meanK, label="K'", linewidth=3) 
        plot!(x_axis, meanB, label="B'", linewidth=3)
        
        plott2 = plot(x_axis, n_share, label="share of exits", linewidth=3)
        
        return ( plot(plott1, plott2, layout=(1,2), size=(1200, 800))   )

    end

    x_size, e_size, _, _ = gridsize()
    plot(xkb_simul(SumPol, Fmat, simn_length = 100000, e_i = e_size-2),
        xkb_simul(SumPol, Fmat, simn_length = 100000, e_i = e_size-1),
        xkb_simul(SumPol, Fmat, simn_length = 100000, e_i = e_size), layout=(4,1), size=(1200, 800))

############ Results: policies ##############
    @elapsed SumPolabl, _, _ = FirmOptim(wage, phi_c = 0)
    @elapsed SumPolcfl, _, _ = FirmOptim(wage, phi_c = 0.7)

    function plotPolicy(SumPolabl, SumPolcfl, e_val)
        
        SPabl = copy(SumPolabl)
        SPcfl = copy(SumPolcfl)
        xvals = string.(Int.(round.(unique(SPabl[:, 1]))))

        n = size(SPabl, 1)
        for s_i in 1:n
        
            SPabl[s_i, 6] + SPabl[s_i, 7] == 1 ? SPabl[s_i, :] .= NaN : SPabl[s_i, :]
            SPcfl[s_i, 6] + SPcfl[s_i, 7] == 1 ? SPcfl[s_i, :] .= NaN : SPcfl[s_i, :]
            
        end
        
        x_size, _, _, _ = gridsize()
        e_ind = e_val*x_size
        
        plott1 = plot(xvals, SPabl[e_ind-x_size+1:e_ind, 3], label="K", linewidth=3, alpha=0.8, xrotation=45)
        plot!(xvals, SPcfl[e_ind-x_size+1:e_ind, 3], label="K0", linewidth=3, linestyle=:dash, alpha=0.8)
        
        plott2 =plot(xvals, SPabl[e_ind-x_size+1:e_ind, 4], label="B", linewidth=3, alpha=0.8, xrotation=45)
        plot!(xvals, SPcfl[e_ind-x_size+1:e_ind, 4], label="B0", linewidth=3, linestyle=:dash, alpha=0.8)
        
        plott3 =plot(xvals, SPabl[e_ind-x_size+1:e_ind, 5], label="X", linewidth=3, alpha=0.8, xrotation=45)
        plot!(xvals, SPcfl[e_ind-x_size+1:e_ind, 5], label="X0", linewidth=3, linestyle=:dash, alpha=0.8)
        
        plott5 =plot(xvals, SPabl[e_ind-x_size+1:e_ind, 9], label="q", linewidth=3, alpha=0.8, xrotation=45)
        plot!(xvals, SPcfl[e_ind-x_size+1:e_ind, 9], label="q0", linewidth=3, linestyle=:dash, alpha=0.8)
        
        plott6 =plot(xvals, SPabl[e_ind-x_size+1:e_ind, 11], label="y", linewidth=3, alpha=0.8, xrotation=45)
        plot!(xvals, SPcfl[e_ind-x_size+1:e_ind, 11], label="y0", linewidth=3, linestyle=:dash, alpha=0.8)
        
        plott7 =plot(xvals, SPabl[e_ind-x_size+1:e_ind, 18]-SPabl[e_ind-x_size+1:e_ind, 13], label="Continuation Value", linewidth=3, alpha=0.8, xrotation=45)
        plot!(xvals, SPcfl[e_ind-x_size+1:e_ind, 18]-SPcfl[e_ind-x_size+1:e_ind, 13], label="Continuation Value0", linewidth=3, linestyle=:dash, alpha=0.8)
        
        return(plot(plott1, plott2, plott3, plott5, plott6, plott7, layout=(3,3), size=(1200, 1100)))
        
                    
    end

    plotPolicy(SumPolabl, SumPolcfl, 11)

############ Results: solving for wage ##############
function FindWage(wage; phi_c)  

    beta = parameters().beta
    SumPol, e_chain, _ = FirmOptim(wage, phi_c = phi_c)

    # productivities set to be equal to the stationary distribution of e_chain    
    e_entry  = reduce(+,stationary_distributions(e_chain))

    # entrant X distribution - x_e = 0  in every case 
    x_vals = unique(SumPol[:, 1])
    zero_index = findall(x -> x == 0.0, x_vals)
    x_entry = zeros(length(x_vals))
    x_entry[zero_index .+ 0] .= 1 # if x = 0 prob = 1 

    # (x,e) are independent, the joint of the two distribution is their product
    xe_entry = [kron(e_entry, x_entry); 0] # also the f0 vector

    # map the entry probabilities to values
    Ve = transpose(xe_entry) * (SumPol[:,end])*beta

    return ( Ve )    

end

#Finding wage given entry cost, using bisection 
tolerance = 1
c_e = 812.0357021040768
# c_e = 3626.022492315388
c_e = 1822.1138130133795


wage_abl = 1
@elapsed SumPol, e_chain, Fmat = FirmOptim(wage_abl, phi_c = 0)
c_e, f0 = EntryValue(SumPol, e_chain) # free entry condition
results_abl = sumSS(SumPol,Fmat,f0)

phi_c = 0.7
@elapsed wage_cfl = find_zero(wage -> FindWage(wage, phi_c = phi_c) - c_e, (0.9, 1.1), Bisection(), rtol=tolerance, verbose=true)
SumPol, e_chain, Fmat = FirmOptim(wage_cfl, phi_c = phi_c)
c_e_cfl, f0 = EntryValue(SumPol, e_chain)
results_cfl = sumSS(SumPol,Fmat,f0)


wage_cfl = 2.052978515625 
wage_cfl = 1.028594970703125  # 1.05841064453125

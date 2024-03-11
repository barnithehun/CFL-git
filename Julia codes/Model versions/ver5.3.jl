########################################################################
###x########### VER 5.3 - trying to functionize ver 5.1 ################
########################################################################

using LinearAlgebra, Statistics, LaTeXStrings, Plots, QuantEcon, Roots
using SparseArrays, ProgressMeter, Dates, XLSX, DataFrames, Distributions


 ####### FIRM OPTIM #######
 function FirmOptim(; 
    # grids sizes - x,k,b should be even numbers!!
        x_size = 60,
        e_size = 7,
        k_size = 40,
        b_size = 40,

    # AR(1) produictivity process
        rho_e = 0.969,
        sigma_e = 0.146,
        nul_e = 1,

    # firm paramaters 
        wage = 1.85,               # will be set by market clearing
        DRS = 0.8,                 # Span of control parameter
        alpha = 1/3*DRS,           # capital share
        nu = 2/3*DRS,              # labour share
        pc = 15,                   # participation cost
        beta = 0.97,               # discount rate
        delta = 0.02,              # depreciation
        pdef_exo = 0.02,           # default probability
        discount = beta*(1-pdef_exo),  # discount parameter
        phi_a = 0.5               # re-sale value of capital
        )

    # functions: contenporaneous optimization
        fn_N(k,e) =  (nu*e*k^alpha / wage)^(1/(1-nu))
        fn_Y(k, e) =  e*k^alpha*fn_N(k,e)^nu
        fn_Pi(k, e) = (1-nu)*fn_Y(k,e)-pc
        fn_X(k,b,e) =  fn_Pi(k, e) + (1-delta) * k - b
        fn_Q(next_k, next_b, pdef) = (next_b == 0) ? 0.97 : (beta * ((1-pdef) * next_b + pdef * min(next_b, phi_a * (1-delta) * next_k))) / next_b
        fn_D(next_k, next_b, x, q) =  x - next_k + q * next_b

    # Setting the state-space
        # productivity process
        e_chain = tauchen(e_size, rho_e, sigma_e, (1-rho_e)*nul_e, 2)
        e_vals = e_chain.state_values
        e_vals = exp.(e_vals)

        # Log-grids
        k_grid = [0;exp.(range(log(10), log(10^5), k_size-1))]
        b_grid = [0;exp.(range(log(10), log(10^5), b_size-1))] # no savings

        x_low = fn_X(k_grid[1],b_grid[end], e_vals[1])
        x_high = fn_X(k_grid[end], b_grid[1], e_vals[end])
        x_grid = [ sort(-exp.(range(log(10), log(-x_low), div(x_size, 2)))[2:end], rev = false); 0; exp.(range(log(10), log(x_high), div(x_size, 2)))] 

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
                    p_trans = e_chain.p[e, next_e_i] 
            
                    # find the second closest
                    if x_next < x_grid[end] && x_next > x_grid[1]
                        
                        if x_next > x_grid[x_close]
                            x_far = x_close + 1
                            else
                            x_far = x_close - 1
                        end
                        
                        x_close_weight = abs(x_next - x_grid[x_far]) / (abs(x_next - x_grid[x_close]) + abs(x_next - x_grid[x_far]))
            
                        # finding the correspoing indicies     
                        xe_close = x_close + (next_e_i-1)*x_size
                        xe_far = x_far + (next_e_i-1)*x_size
                        
                        # filling the transition matrix
                        Q[s_i, a_i, xe_close] = p_trans*x_close_weight
                        Q[s_i, a_i, xe_far] = p_trans*(1-x_close_weight)

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
    q_sa = fill(0.97,n,m)
    pdef_sa = zeros(n,m);
    kbexq_old = zeros(n,4)
    kbexq_new = fill(1,n,4)
    SumPol = zeros(n, 11)
    ################ 
    while !isequal(kbexq_old,kbexq_new)

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

        ddp = DiscreteDP(R, Q, discount);
        results = QuantEcon.solve(ddp, PFI)

        # interpreting results
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

            # Def
            if def == 0
                pdef = pdef_sa[s_i, pol]
                q = q_sa[s_i, pol]
                d = fn_D(k, b, x, q)
            else
                q = 0
                pdef = 0
                d = 0
            end

            # value
            val = values[s_i]
            SumPol[s_i, :] .= [x, e, k, b, next_x, exit, def, pdef, q, d, val]    
        end
        
        ###############################################################################
        # Probability of default and implied q given optimal k', b' in each state
            for s_i in 1:n

            x_i = s_i_vals[s_i, 1]
            e_i = s_i_vals[s_i, 2]

            for a_i in 1:m
                # policies given (x,e)
                next_k = a_vals[a_i,1]
                next_b = a_vals[a_i,2]

                pdef_endo = 0

                for next_e_i in 1:e_size

                    p_trans = e_chain.p[e_i, next_e_i]
                    x_next = fn_X(next_k,next_b,e_vals[next_e_i])

                    x_close = argmin(abs.(x_next .- x_grid))   
                    xe_close = x_close + (next_e_i-1)*x_size
                    next_def_close = a_i_vals[policies[xe_close], 3]

                    if x_next < x_grid[end] && x_next > x_grid[1]
                        
                        if x_next > x_grid[x_close]
                            x_far = x_close + 1
                            else
                            x_far = x_close - 1
                        end
                        
                        x_close_weight = abs(x_next - x_grid[x_far]) / (abs(x_next - x_grid[x_close]) + abs(x_next - x_grid[x_far]))
            
                        # finding the correspoing indicies     
                        xe_far = x_far + (next_e_i-1)*x_size

                        # default policy given (x,e)
                        next_def_far = a_i_vals[policies[xe_far], 3]
                        
                        # update endogeneous default probability
                        pdef_endo += p_trans*(x_close_weight*next_def_close + (1-x_close_weight)*next_def_far)

                    else # x_close_weight = 1
                        pdef_endo += p_trans * next_def_close
                    end  
                
                end 

                # saving results for summary
                pdef = 1 - ((1 - pdef_exo) * (1 - pdef_endo))
                pdef_sa[s_i, a_i] = pdef
                q_sa[s_i,a_i] = fn_Q(next_k, next_b, pdef)

            end
        end

    kbexq_new = SumPol[:, [3, 4, 7, 9]]

    end

    return ( SumPol, e_chain )

end


####### ENTRANTS #######

# mapping the productivity process
function EntryValue(; wage)  

    SumPol, e_chain = FirmOptim(wage = wage)

    # entrant log-productivity distribution
    #   here, set to be equal to the stationary distribution of e
    #   alternatively, set it independently similar to e_chain (they you prly have to scale up the probabilities)
    e_entry  = reduce(+,stationary_distributions(e_chain))

    # entrant X distribution
    x_vals = unique(SumPol[:, 1])
    x_dist = 200*LogNormal(0.4,0.4)-300

    x_entry = zeros(length(x_vals))
    for i = eachindex(x_vals)
        if i == 1
            x_entry[i] = cdf(x_dist, x_vals[i])
        else
            x_entry[i] = cdf(x_dist, x_vals[i]) - cdf(x_dist, x_vals[i-1])
        end
    end

    # (x,e) are independent, the joint of the two distribution is their product
    xe_entry = kron(x_entry, e_entry)

    # map the entry probabilities to values
    Ve = transpose(xe_entry) * (SumPol[1:(end-1),end]-SumPol[1:(end-1),1])
    return(Ve)    

end

tolerance = 1e-1
# Use find_zero with the Newton method and error tolerance
@elapsed result = find_zero(wage -> EntryValue(; wage) - 5500, (1.0, 4.0), Bisection(), rtol=tolerance, verbose=true)


#=
column_names = [:x, :e, :k, :b, :next_x, :exit, :def, :pdef, :q, :d, :val]
SumPol_df = DataFrame(SumPol, column_names)
time_string = "$(Dates.day(now()))_$(Dates.hour(now()))$(Dates.minute(now()))$(Dates.second(now()))"
XLSX.writetable("$time_string.xlsx", sheetname="sheet1", SumPol_df)
#=
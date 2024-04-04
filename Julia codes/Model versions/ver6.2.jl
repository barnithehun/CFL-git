########################################################################
################## VER 6.2 -  just a reparametrization #################
########################################################################

using LinearAlgebra, Statistics, LaTeXStrings, Plots, QuantEcon, Roots, NamedArrays
using SparseArrays, Dates, XLSX, DataFrames, Distributions, Random

function gridsize()
    # grids sizes - x,k,b should be even numbers!!
    return(
        x_size = 60,
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
    pc = 15
    beta = 0.96
    delta = 0.06
    pdef_exo = 0.04
    discount = beta * (1 - pdef_exo)
    phi_a = 0.4
    tauchen_sd = 4
    phi_c = 0.5               # variable cost of reorganization                
    kappa = 0.3               # capital recovery rate of CFL debt
    zeta = 10^4            # fixed cost of reorganization
    tau_vec = 0:1              # vector of CFL reliances - optimal will be one or zero anyways

    return (rho_e = rho_e, sigma_e = sigma_e, nul_e = nul_e, alpha = alpha,
            nu = nu, pc = pc, beta = beta, delta = delta, pdef_exo = pdef_exo,
            discount = discount, phi_a = phi_a, tauchen_sd = tauchen_sd, 
            phi_c = phi_c, kappa = kappa, zeta = zeta, tau_vec = tau_vec)
end


wage = 2

rho_e, sigma_e, nul_e, alpha, nu, pc, beta, delta, pdef_exo, discount, phi_a, tauchen_sd, phi_c, kappa, zeta, tau_vec = parameters()

# calling grid size
x_size, e_size, k_size, b_size = gridsize()           

e_chain = tauchen(e_size, rho_e, sigma_e,0, tauchen_sd)
e_vals = exp.(e_chain.state_values)


# functions: contenporaneous optimization
fn_N(k,e) =  (nu*e*k^alpha / wage)^(1/(1-nu));      
fn_Y(k, e) =  e*k^alpha*fn_N(k,e)^nu;                 
fn_Pi(k, e) = (1-nu)*fn_Y(k,e)-pc;   

# functions: cash on hand and future values
fn_X(k,b,e) =  fn_Pi(k, e) + (1-delta) * k - b;

# Optimal CFL reliance and interest rate
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

# dividends, defined through x 
fn_D(next_k, next_b, x, q) =  x - next_k + q * next_b;

# liquidation decision
fn_chi(k,val) = Int(phi_a*(1-delta)*k >= phi_c*val - zeta)

# Creating grids - should contain 0 in every case 
k_grid = range(0, 10^5, k_size)
b_grid = range(0, 10^5, b_size)
#b_grid =  [range(-10^5, 0, div(b_size, 2))[1:end-1]; range(0, 10^5, div(b_size, 2) +1 )]

x_low = fn_X(k_grid[1],b_grid[end], e_vals[1])
x_high = fn_X(k_grid[end], b_grid[1], e_vals[end])
x_grid =  [range(x_low, 0, div(x_size, 2))[1:end-1]; range(0, x_high, div(x_size, 2) +1 )]


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
@elapsed for a_i in 1:m
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
                    
                    close_weight = abs(x_next - x_grid[x_far]) / (abs(x_next - x_grid[x_close]) + abs(x_next - x_grid[x_far]))
        
                    # finding the correspoing indicies     
                    xe_close = x_close + (next_e_i-1)*x_size
                    xe_far = x_far + (next_e_i-1)*x_size
                    
                    # filling the transition matrix
                    Q[s_i, a_i, xe_close] = p_trans*close_weight
                    Q[s_i, a_i, xe_far] = p_trans*(1-close_weight)
        
                end
        
                # deal with the end points
                if x_next >= x_grid[end] || x_next <= x_grid[1]
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
q_sa = fill(0.97,n,m);
pdef_sa = zeros(n,m);
gam_sa = zeros(n,m);
Pi_liq_sa = zeros(n,m);
Pi_reo_sa =  zeros(n,m);
tau_sa = zeros(n,m);

kbexqt_old = zeros(n,4);
kbexqt_new = similar(kbexqt_old);

nvar_sumres = 15
sumres = zeros(n, nvar_sumres)
################ 
iter = 0
@elapsed while !isequal(kbexqt_old,kbexqt_new)

   iter = iter + 1
   println(iter)
   kbexqt_old = kbexqt_new

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
                R[s_i, a_i] = -1000   # -5000 if you want anyone to quit on its own
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
        
        if def == 0
            q = q_sa[s_i, pol]
            tau = tau_sa[s_i, pol]
            pdef = pdef_sa[s_i, pol]
            gam = gam_sa[s_i, pol]    
            d = fn_D(k, b, x, q)
            Pi_liq = Pi_liq_sa[s_i, pol] 
            Pi_reo = Pi_reo_sa[s_i, pol] 
        else
            q = pdef = d = gam = Pi_liq = Pi_reo = tau = 0
        end

        val = values[s_i]

        sumres[s_i, :] .= [x, e, k, b, next_x, exit, def, pdef, q, gam, d, val, Pi_liq, Pi_reo, tau]
    end


    ###############################################################################
    # Calculates probability of default given optimal k', b' policies that are given by x, e
    # Then it calculates q based on that
    for s_i in 1:n

        x_i = s_i_vals[s_i, 1]
        e_i = s_i_vals[s_i, 2]
    
        for a_i in 1:m
            next_k = a_vals[a_i, 1]
            next_b = a_vals[a_i, 2]
    
            pdef_endo = 0
            gam = 0
            Pi_reo = 0
            Pi_liq = phi_a*(1-delta)*next_k
    
            for next_e_i in 1:e_size
                
                p_trans = e_chain.p[e_i, next_e_i]
                x_next = fn_X(next_k, next_b, e_vals[next_e_i])

                x_index = argmin(abs.(x_next .- x_grid))
                xe_index = x_index + (next_e_i - 1) * x_size
    
                val = values[xe_index]
                pol = policies[xe_index]
                next_def = a_i_vals[pol, 3]
    
                pdef_endo += p_trans * next_def
                gam += p_trans * fn_chi(next_k, val)
                Pi_reo += p_trans * max(phi_c*val - zeta, 0)
            end
    
            pdef = 1 - ((1 - pdef_exo) * (1 - pdef_endo))
            q, tau = fn_Tau_Q(pdef, gam, Pi_liq, Pi_reo, next_b, tau_vec)


            # saving results for summary
            pdef_sa[s_i, a_i] = pdef
            gam_sa[s_i, a_i] = gam
            Pi_liq_sa[s_i, a_i] = Pi_liq
            Pi_reo_sa[s_i, a_i] = Pi_reo
            q_sa[s_i,a_i] = q 
            tau_sa[s_i,a_i] = tau


        end
    end

 kbexqt_new = sumres[:, [3, 4, 6, 15]]

end



##### Print w = 1 results
function PrintPol()    
    column_names = [:x, :e, :k, :b, :next_x, :exit, :def, :pdef, :q, :gam, :d, :val, :Pi_liq, :Pi_reo, :tau]
    SumPol_df = DataFrame(sumres, column_names)
    time_string = "$(Dates.day(now()))_$(Dates.hour(now()))_$(Dates.minute(now()))_$(Dates.second(now()))"
    XLSX.writetable("$time_string.xlsx", sheetname="sheet1", SumPol_df)
end
PrintPol()    

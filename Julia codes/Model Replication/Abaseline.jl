###########################################################################
############### VER 10.2 - Optimal liquidation decision  ##################
###########################################################################

using LinearAlgebra, Statistics, LaTeXStrings, Plots, QuantEcon, Roots, NamedArrays, SparseArrays, Dates, XLSX, DataFrames, Distributions, Random, Optim, Measures

cd(@__DIR__)   # REPLICATION: run from this folder; includes below point to the copies in ./Functions, outputs are saved here

################ Importing Result and Diagnostic Functions ###################
include("Functions/dynsim.jl")
include("Functions/dynsim2.jl")
include("Functions/PrintPolOld.jl")
include("Functions/plotPol.jl")
include("Functions/StatDist_plot.jl")
include("Functions/sumSS.jl")
include("Functions/sumSSsme.jl")
include("Functions/PolperSize.jl")
include("Functions/GamPol.jl")
include("Functions/ZetaGam.jl")
include("Functions/ZetaGamPlot.jl")
include("Functions/CapAlloc.jl")
################ Importing Model Functions ###################
include("Functions/FCmodel/fn_Tau_Q.jl")
include("Functions/FCmodel/EntryValue.jl")
include("Functions/FCmodel/stat_dist.jl")
include("Functions/FCmodel/FindWage.jl")
include("Functions/FCmodel/GridMake.jl")
include("Functions/FCmodel/ErrorFunc.jl")
##################################################################

function gridsize()
    # grid sizes - x, k, b should be even numbers!!
    x_size::Int = 46
    e_size::Int = 23
    k_size::Int = 46
    b_size::Int = 46
    return (x_size = x_size, e_size = e_size, k_size = k_size, b_size = b_size)
end

function parameters()
    rho_e::Float64 = 0.969
    sigma_e::Float64 = 0.146
    nul_e::Int = 1
    DRS::Float64 = 0.75
    alpha::Float64 = 1/3 * DRS
    nu::Float64 = 2/3 * DRS
    pc::Float64 = 42.03135921
    beta::Float64 = 0.96
    delta::Float64 = 0.065
    pdef_exo_l::Float64 =  0.029176179
    pdef_exo_s::Float64 =  pdef_exo_l
    discount::Float64 = beta
    phi_a::Float64 =   0.4
    phi_af::Float64 =   0.4
    tauchen_sd::Float64 = 3

    kappa::Float64 = 0.1       # capital recovery rate of CFL debt
    tau_vec::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}} = 0:1.0:1  # vector of CFL reliances
    zeta_L::Float64 = 0
    phi_c_hh::Float64 = 0.803835125
    Fcut::Float64 = 2000

    return (rho_e = rho_e, sigma_e = sigma_e, nul_e = nul_e, alpha = alpha,
            nu = nu, pc = pc, beta = beta, delta = delta, pdef_exo_s = pdef_exo_s, pdef_exo_l = pdef_exo_l,
            discount = discount, phi_a = phi_a, phi_af = phi_af, tauchen_sd = tauchen_sd,
            kappa = kappa, zeta_L = zeta_L, tau_vec = tau_vec,
             phi_c_hh = phi_c_hh, Fcut = Fcut)
end

####### TRANSITIONS AND POLICY ITERATION #######
# Qv!: out[s,a] = E[ v(s') | s, a ], the same as Q*v in QuantEcon's PFI, computed from Qtab without the dense Q.
#   v is the current value-function guess; it changes every PFI round, while Qtab stays fixed.
#   1. average v over the two x' grid points   2. average over e' with e_ptrans   3. mix in P_x * v[default]
# bellman! then uses vals = R + beta*out to pick the new policy.
function Qv!(out, Qtab, v)
    n, m = size(out)
    x_size, e_size = Qtab.x_size, Qtab.e_size

    # expected v over x' for each action and next productivity e'
    EV = zeros(e_size, m)
    for a_i in 1:m
        Qtab.reg_a[a_i] || continue
        for next_e_i in 1:e_size
            base = (next_e_i - 1) * x_size
            x_close = Qtab.x_close_ea[next_e_i, a_i]
            x_far = Qtab.x_far_ea[next_e_i, a_i]
            close_weight = Qtab.close_weight_ea[next_e_i, a_i]
            EV[next_e_i, a_i] = x_far == 0 ? v[x_close + base] :
                close_weight * v[x_close + base] + (1 - close_weight) * v[x_far + base]
        end
    end

    # expectation over e' given e, and the exogenous default
    EV_e = Qtab.e_ptrans * EV
    for a_i in 1:m, s_i in 1:n
        if s_i == n || !Qtab.reg_a[a_i]
            out[s_i, a_i] = v[n]
        else
            pdef_exo = Qtab.Pdef_a[a_i]
            out[s_i, a_i] = (1 - pdef_exo) * EV_e[Qtab.e_of_s[s_i], a_i] + pdef_exo * v[n]
        end
    end
    return out
end

# Q_sigma_mat builds the transition matrix under the current policy sigma:
# Q_sigma[s, s'] = probability of moving from state s to state s' next period when the firm in state s takes action sigma(s)
# Same as in the original only now the transition matrix is built from the compact Qtab instead of the dense Q.
function Q_sigma_mat(Qtab, sigma)
    n = length(sigma)
    x_size, e_size = Qtab.x_size, Qtab.e_size
    Q_sigma = zeros(n, n)
    for s_i in 1:n
        a_i = sigma[s_i]
        if s_i == n || !Qtab.reg_a[a_i]
            Q_sigma[s_i, n] = 1
        else
            pdef_exo = Qtab.Pdef_a[a_i]
            Q_sigma[s_i, n] = pdef_exo
            for next_e_i in 1:e_size
                base = (next_e_i - 1) * x_size
                p_trans = Qtab.e_ptrans[Qtab.e_of_s[s_i], next_e_i] * (1 - pdef_exo)
                close_weight = Qtab.close_weight_ea[next_e_i, a_i]
                Q_sigma[s_i, Qtab.x_close_ea[next_e_i, a_i] + base] = p_trans * close_weight
                x_far = Qtab.x_far_ea[next_e_i, a_i]
                if x_far != 0
                    Q_sigma[s_i, x_far + base] = p_trans * (1 - close_weight)
                end
            end
        end
    end
    return Q_sigma
end

# bellman!: one policy-improvement step.
#   vals[s,a] = R[s,a] + beta * E[ v(s') | s, a ]   (Qv! fills the expectation in place)
#   sigma[s]  = first a with the highest vals[s,a]   -> the new policy
#   Tv[s]     = that highest value; only the running max for the search, not used afterwards
function bellman!(vals, Tv, sigma, R, Qtab, v, beta)
    Qv!(vals, Qtab, v)
    vals .= R .+ beta .* vals
    n, m = size(vals)
    for s_i in 1:n
        Tv[s_i] = vals[s_i, 1]
        sigma[s_i] = 1
    end
    for a_i in 2:m, s_i in 1:n
        if vals[s_i, a_i] > Tv[s_i]
            Tv[s_i] = vals[s_i, a_i]
            sigma[s_i] = a_i
        end
    end
    return Tv, sigma
end

# PFI_solve: QuantEcon's PFI, using Qtab instead of the dense Q.
#   Evaluate: v = (I - beta*Q_sigma) \ R_sigma, the value of following sigma forever
#   Improve:  sigma = argmax_a R + beta*E[v'], via bellman!/Qv!
#   Repeat until sigma no longer changes.
# Same steps as QuantEcon.solve(DiscreteDP(R, Q, beta), PFI) in QuantEcon 0.19.0
function PFI_solve(R, Qtab, beta; max_iter = 250)
    n, m = size(R)
    vals = similar(R)

    # initial value v(s) = max_a R(s,a) and its greedy policy
    v = vec(maximum(R, dims = 2))
    Tv = similar(v)
    sigma = zeros(Int, n)
    bellman!(vals, Tv, sigma, R, Qtab, v, beta)

    old_sigma = copy(sigma)
    num_iter = 0
    for i in 1:max_iter
        # policy evaluation: value of following sigma forever
        R_sigma = [R[s_i, sigma[s_i]] for s_i in 1:n]
        Q_sigma = Q_sigma_mat(Qtab, sigma)
        v = (I - beta * Q_sigma) \ R_sigma

        # policy improvement: v-greedy policy
        bellman!(vals, Tv, sigma, R, Qtab, v, beta)

        num_iter += 1
        if old_sigma == sigma
            break
        end
        copyto!(old_sigma, sigma)
    end

    return (v = v, sigma = sigma, num_iter = num_iter)
end

####### FIRM OPTIM #######
function FirmOptim(wage, phi_c, zeta_Rl)

    rho_e, sigma_e, nul_e, alpha, nu, pc, beta, delta, pdef_exo_s, pdef_exo_l, discount, phi_a, phi_af, tauchen_sd, kappa,  zeta_L, tau_vec, phi_c_hh, Fcut = parameters()

    # calling grid size
    x_size, e_size, k_size, b_size = gridsize()

    # setting optimization functions
    fn_L(k,e) =  (nu*e*k^alpha / wage)^(1/(1-nu))
    fn_Y(k, e) =  e*k^alpha*fn_L(k,e)^nu
    fn_Pi(k, e) = (1-nu)*fn_Y(k,e)-pc
    fn_X(k,b,e) =  fn_Pi(k, e) + (1-delta) * k - b
    fn_D(next_k, next_b, x, q) =  x - next_k + q * next_b
    
    # Liquidation policies
    fn_Gam(val,k,zeta_R) = Int( (phi_a*(1-delta)*k) >= (phi_c_hh*val - zeta_R)) # optimal liquidation

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
    # Transition probabilities. Q[s,a,:] depends on the state only through e, so instead of the dense
    # Q(n,m,n) array we store where x'(k',b',e') falls on the x-grid, for each action a and next productivity e'.
    # For a regular state (x,e) and a regular action (k',b'):
    #   Q[s, a, xe_close] = e_ptrans[e, e'] * (1-pdef_exo) * close_weight
    #   Q[s, a, xe_far]   = e_ptrans[e, e'] * (1-pdef_exo) * (1-close_weight)
    #   Q[s, a, end]      = pdef_exo   (exogenous default)

    # From the default state, and after the default or exit action, the next state is default with probability 1
    e_of_s = s_i_vals[:, 2]                                    # productivity index of each state
    reg_a = (a_i_vals[:, 3] .== 0) .& (a_i_vals[:, 4] .== 0)   # actions that are neither default nor exit
    Pdef_a = [a_vals[a_i, 1] >= Fcut ? pdef_exo_l : pdef_exo_s for a_i in 1:m]
    x_close_ea = zeros(Int, e_size, m)
    x_far_ea = zeros(Int, e_size, m)      # stays 0 if x' is off the grid, then all weight is on x_close
    close_weight_ea = ones(e_size, m)

    for a_i in 1:m
        reg_a[a_i] || continue
        b = a_vals[a_i, 2]   
        k = a_vals[a_i, 1]   

        for next_e_i in 1:e_size  
            # next period cash on hand x'(b',k',e')
            x_next = fn_X(k,b,e_vals[next_e_i])
            # where x falls on the grid - closer
            x_close = argmin(abs.(x_next .- x_grid))
            x_close_ea[next_e_i, a_i] = x_close

            # find the second closest
            if x_next < x_grid[end] && x_next > x_grid[1]
                x_far = x_next > x_grid[x_close] ? x_close + 1 : x_close - 1            
                close_weight = abs(x_next - x_grid[x_far]) / (abs(x_next - x_grid[x_close]) + abs(x_next - x_grid[x_far]))
                x_far_ea[next_e_i, a_i] = x_far
                close_weight_ea[next_e_i, a_i] = close_weight
            end
        end
    end

    Qtab = (x_size = x_size, e_size = e_size, e_ptrans = e_ptrans, e_of_s = e_of_s, Pdef_a = Pdef_a,
            reg_a = reg_a, x_close_ea = x_close_ea, x_far_ea = x_far_ea, close_weight_ea = close_weight_ea)

    # initial (!) endogenous default probability for each state
    kbexq_old::Array{Float64, 2} = zeros(n, 4)
    kbexq_new::Array{Float64, 2} = fill(1.0, n, 4)
    SumPol::Array{Float64, 2} = zeros(n, 18)
    
    pdef_ea::Array{Float64, 2} = zeros(e_size, m)
    q_ea::Array{Float64, 2} = zeros(e_size, m)
    gam_ea::Array{Float64, 2} = zeros(e_size, m)
    Pi_liq_ea::Array{Float64, 2} = zeros(e_size, m)
    Pi_reo_ea::Array{Float64, 2} = zeros(e_size, m)
    tau_ea::Array{Float64, 2} = zeros(e_size, m)
    iter::Int = 0
    
    ############### 
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
                    q = q_ea[e_of_s[s_i],a_i]
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

        results = PFI_solve(R, Qtab, discount)

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
                e_i = e_of_s[s_i]
                q = q_ea[e_i, pol]
                tau = tau_ea[e_i, pol]
                gam = gam_ea[e_i, pol]
                pdef = pdef_ea[e_i, pol]
                Pi_liq = Pi_liq_ea[e_i, pol]
                Pi_reo = Pi_reo_ea[e_i, pol]
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
        # These depend on the state only through e: they are computed for each (e, a)
        for e_i in 1:e_size
    
            for a_i in 1:m
                # policies given (x,e)
                next_k = a_vals[a_i,1]
                next_b = a_vals[a_i,2]

                # exo probability of default and reorganization costs change with size
                zeta_R = zeta_Rl
                pdef_exo = next_k >= Fcut ? pdef_exo_l : pdef_exo_s
    
                pdef_endo = 0
                gam = 0
                Pi_reo = 0
                Pi_liq = max(0, phi_af*(1-delta)*next_k - zeta_L)

                # ChatGPT made this loop more efficient, not a 100% what it does but it cuts comp time in half - see older versions for the original
                for next_e_i in 1:e_size
    
                    p_trans = e_ptrans[e_i, next_e_i]
                    x_next  = fn_X(next_k, next_b, e_vals[next_e_i])
                    base    = (next_e_i - 1) * x_size

                    j = searchsortedlast(x_grid, x_next)  # j ∈ 0:x_size

                    if 1 ≤ j < x_size
                        t   = (x_next - x_grid[j]) / (x_grid[j+1] - x_grid[j])
                        xe0 = j + base
                        xe1 = j + 1 + base

                        @inbounds begin
                            val      = (1 - t) * Values[xe0] + t * Values[xe1]
                            next_def = (1 - t) * a_i_vals[policies[xe0], 3] +
                                    t * a_i_vals[policies[xe1],  3]
                        end

                    else

                        # clamp to nearest edge
                        i  = (j ≤ 0) ? 1 : x_size
                        xe = i + base
                        
                        @inbounds begin
                            val      = Values[xe]
                            next_def = a_i_vals[policies[xe], 3]
                        end
                    end

                    pdef_endo += p_trans * next_def
                    gam       += p_trans * fn_Gam(val, next_k, zeta_R)
                    Pi_reo    += p_trans * (phi_c * val)
                end

                # now default is a decision and a shock at the same time
                # prob. of either pdef of pdef exo occurs- these are independent events
                pdef = pdef_endo + pdef_exo - pdef_endo * pdef_exo

                # adjusting the Pi_reo to the possibility of exo. default shock
                Pi_reo = Pi_reo * (1-pdef_exo)

                # q and tau
                q, tau = fn_Tau_Q(pdef, gam, Pi_liq, Pi_reo, next_b, tau_vec)

                # saving results for summary
                pdef_ea[e_i, a_i] = pdef
                q_ea[e_i,a_i] = q
                tau_ea[e_i,a_i] = tau
                gam_ea[e_i, a_i] = gam
                Pi_liq_ea[e_i, a_i] = Pi_liq
                Pi_reo_ea[e_i, a_i] = Pi_reo
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

        # probability that you will end up in default state
        pdef_exo = next_k >= Fcut ? pdef_exo_l : pdef_exo_s
        Fmat_bottom[s_i] = pdef_exo + SumPol[s_i, 8] * SumPol[s_i, 14]

        e_i = Int(floor( (s_i-1) / x_size) + 1) 
        for next_e_i in 1:e_size

            p_trans = e_ptrans[e_i, next_e_i] * (1-(pdef_exo + SumPol[s_i, 8] * SumPol[s_i, 14]))
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

############ Results: ABL vs CFL  ##############
# calibration values
wage = 1;
phi_c =  0.360542276
zeta_Rl =  2753.229657
@elapsed SumPol, e_chain, Fmat = FirmOptim(wage, phi_c, zeta_Rl)
c_e, f0 = EntryValue(SumPol, e_chain) ;
mu, m, xpol = stat_dist(SumPol, Fmat, f0);

zeta_Rl = 1350
@elapsed SumPol0, e_chain0, Fmat0 = FirmOptim(1.0175, phi_c, zeta_Rl)  # perfcred w is about 1.14
c_e0, f00 = EntryValue(SumPol0, e_chain0); 
mu0, m0, xpol0 = stat_dist(SumPol0, Fmat0, f00);


############ Core results Results: stationary distributions ############
println("The relative productivity of the ABL case is: ", round(c_e0 / c_e , digits=3))

sumSSsme(SumPol, Fmat,f0)
sumSSsme(SumPol0, Fmat0,f00)

sumSS(SumPol,Fmat,f0)
sumSS(SumPol0,Fmat0,f00)


############ Untargeted Moments ############
binnum = 10
# Pdef and CFL reliance
Ushape(binnum, SumPol, mu) # in calculating tau, this contains firms with Pdef = 0
# Gamma and CFL reliance
Xcross(binnum, SumPol, mu) # in calculating tau, this does not contain firms with Pdef = 0
# Average Debt policy and Interest rate across firm sizes
QBplot(binnum, SumPol, SumPol0, mu, mu0) 
# Capital Allocation 
CapAlloc(SumPol, SumPol0, mu, mu0) 


# Q - against Leverage this plot works well with x_size == 46 and e_size == 23 
p1, p2 = DebtScedule(SumPol, 25, 17; phi_c = phi_c)
p3, p4 = DebtScedule(SumPol, 25, 20; phi_c = phi_c)
p5, p6 = DebtScedule(SumPol, 25, 22; phi_c = phi_c)
plot(p3, p4, p5, p6, layout = (2, 2), size = (900, 700))
#savefig(plot(p3, p4, p5, p6, layout = (2, 2), size = (900, 700) ), "DebtScedule.png")

############ Results: dynamics simulations ##############
_, e_size, _, _ = gridsize()
# shows the average policies of firms of a given productivity across their lifecycle
dynsim2(e_size - 6, 10000)
dynsim2(e_size - 3, 10000)
dynsim2(e_size - 1, 10000)
# savefig(dynsim2(e_size - 2, 100000), "Simul.png")

# shows optimal policies of firms with a certain productivity, across x-sizes
plotPol(SumPol0, SumPol, e_size-6)

############ Firm states and policies in stationary equilibrium ###########
binnum = 20
plot(plotPDF(binnum, 'k', SumPol), plotPDF(binnum, 'b', SumPol), plotPDF(binnum, 'l', SumPol),
    plotPDF(binnum, 'y', SumPol), plotPDF(binnum, 'p', SumPol), plotPDF(binnum, 'v', SumPol), layout=(2,3), size=(1200, 800))
plot(plotCDF(binnum, 'k', SumPol), plotCDF(binnum, 'b', SumPol), plotCDF(binnum, 'l', SumPol),
    plotCDF(binnum, 'y', SumPol), plotCDF(binnum, 'p', SumPol), plotCDF(binnum, 'v', SumPol), layout=(2,3), size=(1200, 800))

plotXE(SumPol, mu, e_chain)     
plotTauDist(SumPol)


# Q - against Gam this plot works well with x_size == 44 and e_size == 27 
pGamPol = plot(GamPol(SumPol, 15, 21), GamPol(SumPol, 26, 21), GamPol(SumPol, 35, 21), layout = (1, 3), size = (900, 300) , margin = 5mm)
# Plots.savefig(pGamPol, "GamPol.png")

# Reform effects
include("Functions/improve.jl")
reformplot = plotPolHeatmaps(SumPol, SumPol0, mu, 1)
# Plots.savefig(reformplot, "large_reform.png")
# Plots.savefig(reformplot, "reform.png")


### Optimal policies over different reorganization costs ###
zeta_R_vec = 0:1000:8000
# SumPolZeta = ZetaGam(zeta_R_vec)  # runtime around 1 hour 
ZetaGamPlot(SumPolZeta, zeta_R_vec, 23, 15)

###### GENERAL EQUILIBRIUM: Finding wage given entry cost, using bisection ####
tolerance = 1
wage_0 = 1
phi_c =  0.360542276
zeta_Rl =  2753.229657

SumPol, e_chain, Fmat = FirmOptim(wage_0, phi_c, zeta_Rl)
c_e, f0 = EntryValue(SumPol, e_chain)
results_baseline = sumSS(SumPol,Fmat,f0)
c_e_baseline = c_e

#zeta_Rl = 1350
@elapsed wage_alt = find_zero(wage -> FindWage(wage, phi_c, zeta_Rl) - c_e_baseline, (1.1, 1.3), Bisection(), rtol=tolerance, verbose=true)

#  wage_alt = 1.017549201334813 - where zeta_Rl  = zeta_Rl/2
#  wage_alt =  - where zeta_Rl  = zeta_Rl/10

#  wage_inf =  - where zeta_Rl  = inf
#  wage_alt =  - perfect credit economy 
wage_alt

SumPol0, e_chain0, Fmat0 = FirmOptim(wage_alt, phi_c, zeta_Rl)
c_e0, f00 = EntryValue(SumPol0, e_chain0)

sumSS(SumPol0,Fmat0,f00) # wage_alt = 1.025
sumSSsme(SumPol0,Fmat0,f00)

########## CALIBRATION ################
include("Functions/FCmodel/ErrorFunc.jl")
include("Functions/FCmodel/FirmOptim_Ext.jl")

# initital values ipc, ipdef_exo_l, ipdef_exo_s, izeta_Rl, izeta_Rs, iphi_a, iphi_c, iFcut)
initvec = [ 28.56, 0.02343, 2759.45, 0.75978, 0.393]
options = Optim.Options(
    iterations = 200,
    time_limit = 8*3600.0,
    f_tol = 0.002,
    show_trace = true,
    store_trace = true,
    show_every = 1)
@elapsed CalRes = optimize(ErrorFunc1, initvec, NelderMead(), options)

println(CalRes)


# CFshare_calib = CalRes.minimizer
CFrel_calib = CalRes.minimizer

ErrorFunc2(initvec)
ErrorFunc2(CalRes.minimizer)

# Wrap CFrel_calib into a DataFrame
CFrel_df = DataFrame(CFcal = CFrel_calib)
# Timestamped filename
time_string = "ESTPARAM_$(Dates.day(now()))_$(Dates.hour(now()))_$(Dates.minute(now()))"
# Write to Excel
XLSX.writetable("$time_string.xlsx", sheetname = "sheet1", CFrel_df)



###### GENERAL EQUILIBRIUM: Shifting Costs ####
tolerance = 1
wage_0 = 1
phi_c = 0.2810
zeta_Rl = 2590

SumPol, e_chain, Fmat = FirmOptim(wage_0, phi_c, zeta_Rl)
c_e, f0 = EntryValue(SumPol, e_chain)
results_baseline = sumSS(SumPol,Fmat,f0)
c_e_baseline = c_e

zeta_Rl = 1295
phi_c = 0.14
@elapsed wage_alt = find_zero(wage -> FindWage(wage, phi_c, zeta_Rl) - c_e_baseline, (0.99, 1.02), Bisection(), rtol=tolerance, verbose=true)


SumPol0, e_chain0, Fmat0 = FirmOptim(wage_alt, phi_c, zeta_Rl)
c_e0, f00 = EntryValue(SumPol0, e_chain0)

sumSS(SumPol0,Fmat0,f00) # wage_alt = 1.025
sumSSsme(SumPol0,Fmat0,f00)

sumSS(SumPol,Fmat,f0) # wage = 1
sumSS(SumPol0,Fmat0,f00) # wage_alt = 1.025	


####### Finding zeta_L for a given target liquidation probability ####
include("Functions/FindZeta.jl")

tolerance = 1
wage = 1
phi_c =  0.360542276
zeta_Rl =  2753.229657

FindZeta(wage, phi_c, zeta_Rl)
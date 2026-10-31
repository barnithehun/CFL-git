###########################################################################
############### VER 10.2 - Optimal liquidation decision  ##################
###########################################################################

using LinearAlgebra, Statistics, LaTeXStrings, Plots, QuantEcon, Roots, NamedArrays, SparseArrays, Dates, XLSX, DataFrames, Distributions, Random, Optim, Measures, CSV

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
include("Functions/FCmodel/Calibrate.jl")
##################################################################

function gridsize()
    # grid sizes - x, k, b should be even numbers!!
    x_size::Int = 46
    e_size::Int = 23
    k_size::Int = 46
    b_size::Int = 46
    return (x_size = x_size, e_size = e_size, k_size = k_size, b_size = b_size)
end

# Overrides: parameters() returns these values in place of the defaults written below.
# Only with_params, set_params! and reset_params! (Functions/FCmodel/Calibrate.jl) change it; when it is empty, parameters() is unchanged.
const PARAM_OVERRIDE = Dict{Symbol, Float64}()
par(name, default) = get(PARAM_OVERRIDE, name, default)

function parameters()
    rho_e::Float64 = par(:rho_e, 0.969)
    sigma_e::Float64 = par(:sigma_e, 0.146)
    nul_e::Int = par(:nul_e, 1)
    DRS::Float64 = par(:DRS, 0.75)
    alpha::Float64 = par(:alpha, 1/3 * DRS)
    nu::Float64 = par(:nu, 2/3 * DRS)
    pc::Float64 = par(:pc, 42.03135921)
    beta::Float64 = par(:beta, 0.96)
    delta::Float64 = par(:delta, 0.065)
    pdef_exo_l::Float64 = par(:pdef_exo_l, 0.029176179)
    pdef_exo_s::Float64 = par(:pdef_exo_s, pdef_exo_l)
    discount::Float64 = par(:discount, beta)
    phi_a::Float64 = par(:phi_a, 0.4)
    phi_af::Float64 = par(:phi_af, 0.4)
    tauchen_sd::Float64 = par(:tauchen_sd, 3)

    kappa::Float64 = par(:kappa, 0.1)       # capital recovery rate of CFL debt
    tau_vec::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}} = 0:1.0:1  # vector of CFL reliances
    zeta_L::Float64 = par(:zeta_L, 0)
    phi_c_hh::Float64 = par(:phi_c_hh, 0.803835125)
    Fcut::Float64 = par(:Fcut, 2000)

    # baseline values of the two FirmOptim arguments (phi_c is theta in main5)
    phi_c::Float64 = par(:phi_c, 0.360542276)
    zeta_Rl::Float64 = par(:zeta_Rl, 2753.229657)

    # phi_c, zeta_Rl and DRS come last: FirmOptim, GamPol and CapAlloc2 unpack the first 19 fields by position
    return (rho_e = rho_e, sigma_e = sigma_e, nul_e = nul_e, alpha = alpha,
            nu = nu, pc = pc, beta = beta, delta = delta, pdef_exo_s = pdef_exo_s, pdef_exo_l = pdef_exo_l,
            discount = discount, phi_a = phi_a, phi_af = phi_af, tauchen_sd = tauchen_sd,
            kappa = kappa, zeta_L = zeta_L, tau_vec = tau_vec,
            phi_c_hh = phi_c_hh, Fcut = Fcut, phi_c = phi_c, zeta_Rl = zeta_Rl, DRS = DRS)
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
# AI updateded it for sparse matrix construction, which is more efficient for large state spaces.
function Q_sigma_mat(Qtab, sigma)
    n = length(sigma)
    x_size, e_size = Qtab.x_size, Qtab.e_size
    rows = Int[]; cols = Int[]; probs = Float64[]    # one entry per nonzero Q_sigma[s, s']
    for s_i in 1:n
        a_i = sigma[s_i]
        if s_i == n || !Qtab.reg_a[a_i]
            push!(rows, s_i); push!(cols, n); push!(probs, 1.0)
        else
            pdef_exo = Qtab.Pdef_a[a_i]
            push!(rows, s_i); push!(cols, n); push!(probs, pdef_exo)
            for next_e_i in 1:e_size
                base = (next_e_i - 1) * x_size
                p_trans = Qtab.e_ptrans[Qtab.e_of_s[s_i], next_e_i] * (1 - pdef_exo)
                close_weight = Qtab.close_weight_ea[next_e_i, a_i]
                push!(rows, s_i); push!(cols, Qtab.x_close_ea[next_e_i, a_i] + base); push!(probs, p_trans * close_weight)
                x_far = Qtab.x_far_ea[next_e_i, a_i]
                if x_far != 0
                    push!(rows, s_i); push!(cols, x_far + base); push!(probs, p_trans * (1 - close_weight))
                end
            end
        end
    end
    return sparse(rows, cols, probs, n, n)
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
#   v_init: starting values. FirmOptim passes the previous q iteration's values (warm start); the default is QuantEcon's v(s) = max_a R(s,a).
# Same steps as QuantEcon.solve(DiscreteDP(R, Q, beta), PFI) in QuantEcon 0.19.0
function PFI_solve(R, Qtab, beta; v_init = nothing, max_iter = 250)
    n, m = size(R)
    vals = similar(R)

    # initial value (v_init, or v(s) = max_a R(s,a)) and its greedy policy
    v = v_init === nothing ? vec(maximum(R, dims = 2)) : v_init
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
# phi_c and zeta_Rl default to their values in parameters(); pass them only to solve at other values (e.g. a reform)
function FirmOptim(wage, phi_c = parameters().phi_c, zeta_Rl = parameters().zeta_Rl)

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
    Pdef_a = [next_k >= Fcut ? pdef_exo_l : pdef_exo_s for next_k in a_vals[:, 1]]   # loop over the column, so the comprehension does not capture a_vals (which would box it)
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

    # Tables for the q update: where x'(k',b',e') falls on the x-grid does not change across q iterations
    x_lo_ea = zeros(Int, e_size, m)      # state index of the grid point at or below x' (or of the nearest edge)
    t_ea = zeros(e_size, m)              # interpolation weight on the next grid point
    interp_ea = falses(e_size, m)        # true if x' is inside the grid
    for a_i in 1:m, next_e_i in 1:e_size
        x_next = fn_X(a_vals[a_i, 1], a_vals[a_i, 2], e_vals[next_e_i])
        base = (next_e_i - 1) * x_size
        j = searchsortedlast(x_grid, x_next)  # j ∈ 0:x_size
        if 1 ≤ j < x_size
            t_ea[next_e_i, a_i] = (x_next - x_grid[j]) / (x_grid[j+1] - x_grid[j])
            x_lo_ea[next_e_i, a_i] = j + base
            interp_ea[next_e_i, a_i] = true
        else
            x_lo_ea[next_e_i, a_i] = ((j ≤ 0) ? 1 : x_size) + base   # clamp to nearest edge
        end
    end

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
    Values_prev = nothing   # values from the previous q iteration, the starting point for PFI_solve
    q_tol = 1e-10           # stop when k', b', default are unchanged and q moved by at most q_tol
    
    ############### 
    while !(kbexq_old[:, 1:3] == kbexq_new[:, 1:3] &&
            maximum(abs.(kbexq_new[:, 4] .- kbexq_old[:, 4])) <= q_tol)

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

        results = PFI_solve(R, Qtab, discount; v_init = Values_prev)

        Values = results.v;
        policies = results.sigma;  # optimal policy     
        Values_prev = Values

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

                # Expectation over e' only: given (k', b', e'), x' is known, so only the two x-grid points around x' matter (interpolated)
                for next_e_i in 1:e_size
    
                    p_trans = e_ptrans[e_i, next_e_i]
                    xe0     = x_lo_ea[next_e_i, a_i]   # from the tables built before the q loop

                    if interp_ea[next_e_i, a_i]
                        t   = t_ea[next_e_i, a_i]
                        xe1 = xe0 + 1

                        @inbounds begin
                            val      = (1 - t) * Values[xe0] + t * Values[xe1]
                            next_def = (1 - t) * a_i_vals[policies[xe0], 3] +
                                    t * a_i_vals[policies[xe1],  3]
                        end

                    else

                        # x' off the grid: nearest edge
                        @inbounds begin
                            val      = Values[xe0]
                            next_def = a_i_vals[policies[xe0], 3]
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
                q, tau = fn_Tau_Q(pdef, gam, Pi_liq, Pi_reo, next_b, tau_vec, beta, kappa)

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

########## CALIBRATION ################
# Data targets: :R reads the last table of Data Replication/calib_targets.R (run it first to update calib_targets.csv);
# :manual uses manual_targets below (main5's values). "Share of small firms" sets Fcut after the fit.
targets_source = :R
manual_targets = DataFrame(
    moment   = ["Debt to Collateral", "Interest rate", "Liquidation prob.", "CF-reliance", "CF-based to total debt", "Share of small firms"],
    weighted = [0.52,                 5.30,            0.58,                0.51,          0.82,                     0.66],
    small    = [0.44,                 5.89,            0.84,                0.51,          0.34,                     missing],
    large    = [0.59,                 4.63,            0.20,                0.53,          0.83,                     missing])
targets_source in (:R, :manual) || error("targets_source must be :R or :manual, not $(repr(targets_source))")
data_targets = targets_source == :R ? read_targets(joinpath("..", "..", "Data Replication", "calib_targets.csv")) : manual_targets

# calib_moments: (moment as in calib_targets.R, row in sumSS, row in sumSSsme, error form)
#   :rel = |model - data| / data;  :rel_1m = |model - data| / (1 - data), as ErrorFunc1's AB share;  :none = not targeted
calib_moments = [
    ("Debt to Collateral",     "avg_b2a",     "b2a",         :rel),
    ("Interest rate",          "avg intrate", "intrate",     :rel),
    ("Liquidation prob.",      "avg_liqprob", "liqprob",     :rel),
    ("CF-reliance",            "CF reliance", "CF reliance", :rel),
    ("CF-based to total debt", "CF share",    "CF share",    :rel_1m),
]

# calib_params: name = (start, lower, upper); each trial must lie strictly inside the bounds.
#  Starts are the current values in parameters(). Add a line to calibrate a parameter, delete one to keep its parameters() value.
#  Fcut is not listed: it only splits small and large firms, so it is set after the fit to match the share of small firms.
calib_params = (
    pc         = (parameters().pc,         0.0, Inf),   # fixed cost of operation, c
    pdef_exo_l = (parameters().pdef_exo_l, 0.0, 1.0),   # exogenous distress probability, P_x (pdef_exo_s follows it)
    zeta_Rl    = (parameters().zeta_Rl,    0.0, Inf),   # fixed cost of reorganization, zeta
    phi_c_hh   = (parameters().phi_c_hh,   0.0, 1.0),   # 1 - variable cost of reorganization (main5's phi_c = 1 - phi_c_hh)
    phi_c      = (parameters().phi_c,      0.0, 1.0),   # lender's payoff under reorganization, theta
)

# making that what the lender and the HH gets after reorganization does not add up to more that 1
valid = p -> (1 - p.phi_c_hh) + p.phi_c <= 1

options = Optim.Options(
    iterations = 200,
    time_limit = 3*3600.0,
    f_reltol = 0.002,      # f_tol before Optim 2
    show_trace = true,
    store_trace = true,
    show_every = 1)
# thorough = true runs calibrate_thoroughly instead (see Calibrate.jl): it restarts Nelder-Mead until a restart gains less
#  than 1e-6 (at most 5 restarts), from the start in calib_params and two variants ±variant_pct away, and keeps the best.
#  It takes several times as long as calibrate. starts = [CalRes.values] in its call polishes an earlier result instead.
thorough = false
# saves to calibration_outputs/CALIB_<date>.csv (thorough: CALIB_<date>_thorough.csv and CALIB_<date>_thorough_runs.csv)
@elapsed CalRes = thorough ?
    calibrate_thoroughly(calib_params, calib_moments, data_targets; options = options, valid = valid, variant_pct = 0.10) :
    calibrate(calib_params, calib_moments, data_targets; options = options, valid = valid)

# To use the calibrated values (Fcut included) without pasting them: set_params!(CalRes.values).
# reset_params!() goes back to the values written in parameters().


########## MODEL FIT ################
# Fit at the calibrated values. Without a calibration run, model_fit(NamedTuple(), calib_moments, data_targets)
# gives the fit at the values in parameters().
Fit = model_fit(CalRes.values, calib_moments, data_targets)
print_fit(Fit)



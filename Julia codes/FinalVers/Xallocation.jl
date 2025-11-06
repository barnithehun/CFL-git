###########################################################################
############### VER 10.2 - Optimal liquidation decision  ##################
###########################################################################

using LinearAlgebra, Statistics, LaTeXStrings, Plots, QuantEcon, Roots, NamedArrays, SparseArrays, Dates, XLSX, DataFrames, Distributions, Random, Optim, Measures

################ Importing Result and Diagnostic Functions ###################
include("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Julia codes/Functions/dynsim.jl")
include("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Julia codes/Functions/dynsim2.jl")
include("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Julia codes/Functions/PrintPolOld.jl")
include("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Julia codes/Functions/plotPol.jl")
include("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Julia codes/Functions/StatDist_plot.jl")
include("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Julia codes/Functions/sumSS.jl")
include("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Julia codes/Functions/sumSSsme.jl")
include("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Julia codes/Functions/PolperSize.jl")
include("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Julia codes/Functions/GamPol.jl")
include("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Julia codes/Functions/ZetaGam.jl")
include("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Julia codes/Functions/ZetaGamPlot.jl")
include("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Julia codes/Functions/CapAlloc.jl")
################ Importing Model Functions ###################
include("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Julia codes/Functions/FCmodel/fn_Tau_Q.jl")
include("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Julia codes/Functions/FCmodel/EntryValue.jl")
include("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Julia codes/Functions/FCmodel/stat_dist.jl")
include("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Julia codes/Functions/FCmodel/FindWage.jl")
include("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Julia codes/Functions/FCmodel/GridMake.jl")
include("C:/Users/szjud/OneDrive/Asztali gép/EBCs/CFL-git/Julia codes/Functions/FCmodel/ErrorFunc.jl")
##################################################################

function gridsize()
    # grid sizes - x, k, b should be even numbers!!
    x_size::Int = 56
    e_size::Int = 21
    k_size::Int = 42
    b_size::Int = 42
    return (x_size = x_size, e_size = e_size, k_size = k_size, b_size = b_size)
end

function parameters()
    rho_e::Float64 = 0.969
    sigma_e::Float64 = 0.146
    nul_e::Int = 1
    DRS::Float64 = 0.75
    alpha::Float64 = 1/3 * DRS
    nu::Float64 = 2/3 * DRS
    pc::Float64 = 45.04042
    beta::Float64 = 0.96
    delta::Float64 = 0.065
    pdef_exo_l::Float64 =    0.03161
    pdef_exo_s::Float64 =  pdef_exo_l
    discount::Float64 = beta            
    phi_a::Float64 =   0.4
    phi_af::Float64 =   0.4
    tauchen_sd::Float64 = 3

    kappa::Float64 = 0.1       # capital recovery rate of CFL debt
    tau_vec::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}} = 0:1.0:1  # vector of CFL reliances
    zeta_L::Float64 = 0
    phi_c_hh::Float64 =   0.91508
    Fcut::Float64 = 2000

    return (rho_e = rho_e, sigma_e = sigma_e, nul_e = nul_e, alpha = alpha,
            nu = nu, pc = pc, beta = beta, delta = delta, pdef_exo_s = pdef_exo_s, pdef_exo_l = pdef_exo_l,
            discount = discount, phi_a = phi_a, phi_af = phi_af, tauchen_sd = tauchen_sd,
            kappa = kappa, zeta_L = zeta_L, tau_vec = tau_vec,
             phi_c_hh = phi_c_hh, Fcut = Fcut)
end

####### FIRM OPTIM #######
function FirmOptim(wage, phi_c, zeta_Rl, PerfCred)

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
    SumPol::Array{Float64, 2} = zeros(n, 18)
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

                # if Perfect Credit is turned on negative dividens have no costs (free equity injections)
                if PerfCred == 1 || (PerfCred == 0 && d > 0)
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

                # updating Gamma: probability of liquidation either through endogenous or exogenous
                # gam = (gam*pdef_exo + pdef_endo - ( gam*pdef_exo * pdef_endo)) / pdef

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




# ------------ Setup ------------
# fixed params
wage_baseline = 1.0
phi_c = 0.3340487197
zeta_Rl=3071.212

# scenarios 1–4
scenarios = [
    (id=1, wage=wage_baseline, zeta_Rl=3071.212, perfect=0, label="baseline (ζ)"),
    (id=2, wage=1.020,         zeta_Rl=1500.0,            perfect=0, label="reform (ζ/2)"),
    (id=3, wage=1.020,         zeta_Rl=300.0,             perfect=0, label="reform (ζ/10)"),
    (id=4, wage=1.020,         zeta_Rl=1500.0,            perfect=1, label="perfect credit")
]

# model sizes and parameters
x_size, e_size, k_size, b_size = gridsize()
rho_e, sigma_e, nul_e, alpha, nu, pc, beta, delta, pdef_exo_s, pdef_exo_l, discount, phi_a, phi_af, tauchen_sd, kappa, zeta_L, tau_vec, phi_c_hh, Fcut = parameters()

# productivity grid and "optimal" allocation objects (common across scenarios)
prod_toskip = 0
toskip = x_size * prod_toskip + 1

e_chain = tauchen(e_size, rho_e, sigma_e, (1 - rho_e) * nul_e, tauchen_sd)
e_vals  = exp.(e_chain.state_values)
x       = e_vals[prod_toskip+1:end]

e_stat_full = sum(stationary_distributions(e_chain))
e_stat      = e_stat_full[prod_toskip+1:end]
e_stat     /= sum(e_stat)

k_share_opt = e_vals[prod_toskip+1:end].^(1/(1 - alpha - nu))
k_share_opt ./= sum(k_share_opt)
k_share_opt = e_stat .* k_share_opt
k_share_opt ./= sum(k_share_opt)

# ------------ Run scenarios and collect results ------------
results = Dict{Int, Dict{Symbol,Any}}()

for sc in scenarios
    @elapsed begin
        SumPol, e_chain_s, Fmat = FirmOptim(sc.wage, phi_c, sc.zeta_Rl, sc.perfect)
        c_e, f0  = EntryValue(SumPol, e_chain_s)
        mu, m, xpol = stat_dist(SumPol, Fmat, f0)
        results[sc.id] = Dict(
            :label => sc.label,
            :SumPol => SumPol,
            :mu => mu,
            :m => m,
            :xpol => xpol
        )
    end
end

# ------------ Per-scenario derived objects ------------
for (id, R) in results
    mu = R[:mu]
    Span = toskip:(length(mu) - 1)

    mu_stat = sum(reshape(mu[Span], x_size, :), dims=1)[:]               # firm mass by productivity
    M = sum(mu)                                                           # total mass

    k_num = results[id][:SumPol][Span, 3] .* mu[Span]
    k_share_obs = k_num ./ sum(k_num)
    k_share_obs = sum(reshape(k_share_obs, x_size, :), dims=1)[:]         # by productivity

    R[:mu_stat] = mu_stat
    R[:M] = M
    R[:k_share_obs] = k_share_obs
    R[:check_mu_share_sum] = sum(mu_stat) / M
end

# ------------ Plot 1: Capital allocation (observed vs optimal) ------------
# --- plt1: Capital allocation (skip id=4). id=2 dotted ---
R1 = results[1]
R2 = results[2]
R3 = results[3]

plt1 = plot(x, k_share_opt, label="perfect credit", lw=3,  color=:black, size=(900,600), margin = 5mm)
plot!(plt1, x, R1[:k_share_obs], label=R1[:label], lw=3, color=:blue3)
plot!(plt1, x, R2[:k_share_obs], label=R2[:label], lw=3, color=:red, linestyle=:dashdot, dashes=(1.0, 4.0))
plot!(plt1, x, R3[:k_share_obs], label=R3[:label], linestyle=:dash, lw=3, color=:green)  # dark violet

xlabel!(plt1, "Productivity")
ylabel!(plt1, "Capital Share")
#title!(plt1, "Firm Mass Across Productivity States")
    plot!(plt1, xtickfont=font(10), ytickfont=font(10), guidefont=font(12), legendfont=font(10))

#savefig(plt1, "capital_allocation_4runs.png")

# ------------ Plot 2: Firm distribution by productivity ------------
# --- plt2: Firm distribution (include id=4) ---
R1 = results[1]
R2 = results[2]
R3 = results[3]
R4 = results[4]

plt2 = plot()

plot!(plt2, x, R1[:mu_stat] / sum(R1[:mu_stat]), label=R1[:label], lw=2.5, color=:blue3)
plot!(plt2, x, R2[:mu_stat] / sum(R2[:mu_stat]), label=R2[:label], lw=2.5, color=:red, linestyle=:dashdot, dashes=(1.0, 4.0))
plot!(plt2, x, R3[:mu_stat] / sum(R3[:mu_stat]), label=R3[:label], lw=2.5, color=:green, linestyle=:dash)
plot!(plt2, x, R4[:mu_stat] / sum(R4[:mu_stat]), label=R4[:label], lw=2.5, color=:black)  # violet for scenario 4

xlabel!(plt2, "Productivity")
ylabel!(plt2, "Share of Firm Mass")
#title!(plt2, "Firm Mass Across Productivity States")
plot!(plt2, xtickfont=font(10), ytickfont=font(10),
      guidefont=font(12), legendfont=font(10))

# savefig(plt2, "firm_distributionshare_4runs.png")

# ------------ Plot 3: Firm distribution by productivity ------------
# --- plt2: Firm distribution (include id=4) ---
R1 = results[1]
R2 = results[2]
R3 = results[3]
R4 = results[4]

plt3 = plot()

plot!(plt3, x, R1[:mu_stat] , label=R1[:label], lw=2.5, color=:blue3)
plot!(plt3, x, R2[:mu_stat], label=R2[:label], lw=2.5, color=:red, linestyle=:dashdot, dashes=(1.0, 4.0))
plot!(plt3, x, R3[:mu_stat], label=R3[:label], lw=2.5, color=:green, linestyle=:dash)
plot!(plt3, x, R4[:mu_stat], label=R4[:label], lw=2.5, color=:black)  # violet for scenario 4

xlabel!(plt3, "Productivity")
ylabel!(plt3, "Firm Mass")
#title!(plt3, "Firm Mass Across Productivity States")
plot!(plt3, xtickfont=font(10), ytickfont=font(10),
      guidefont=font(12), legendfont=font(10))

savefig(plt3, "firm_distribution_4runs.png")

plot(plt3, plt1, layout=(1,2), size=(1200, 500), margin=10mm)
savefig("Combined4runs2.png")

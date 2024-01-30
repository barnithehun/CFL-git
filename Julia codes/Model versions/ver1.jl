########################################################################
############ VER 1 - PE model - Fixed Qs - 3 state variables ###########
########################################################################

 using LinearAlgebra, Statistics,  LaTeXStrings, Plots, QuantEcon, SparseArrays

###  State variables (k,b,ε)
###  Action variables (k', b')
###  Only, ε is stochastic, k' and b' are defined fully from previous actions


# Later, put all this into a function, 
# function Firms(; ....)


# grids size  - around 20x20x2 the biggest I can have rn
k_size = 3;
b_size = 3;
# Define idiosyncratic shocks
e_chain = MarkovChain([0.9 0.1; 0.1 0.9], [8; 12])
# two results to extract are, z_chain.p; z_chain.state_values
e_size = length(e_chain.state_values)


# paramaters 
wage = 1.85;            # wage (exo. in the PE model)
alpha = 1/3*0.88;       # capital share
nu = 2/3*0.88;          # labour share
pc = 1.5;               # participation cost
beta = 0.989;           # discount rate
delta = 0.015;          # depreciation
pdef = 0.015;           # default probability
discount = beta*(1-pdef);   # discount parameter
q = 0.94;               # external finance premium 
# having the q like this implies that the firm will probably fully wants to finance from debt
# it is always cheaper than internal finance

# Creating grids
# LIN GRIDS
k_grid = range(3.8*10^4, 10^5, k_size)
b_grid = 0.4*range(-(10^5), 10^5, b_size) 
# note: you need to find a grid where at least one action is feasible in each state

#= LOG GRIDS
 # k-grid
   k_grid = 10 .^ range(0, 6, k_size)
 # this is a log grid between -5*10e5 and 5*10e5, but it is centered around 0
 b_grid = 5 .* [sort(-10 .^ range(0, 5, div(b_size, 2))[2:end], rev = false); 10 .^ range(0, 5, div(b_size, 2) +1 )] 
=#

# functions given (k,ε)
fn_N(k,e) =  (nu*e*k^alpha / wage)^(1/(1-nu));      
fn_Y(k, e) =  e*k^alpha*fn_N(k,e)^nu;                 
fn_Pi(k, e) = (1-nu)*fn_Y(k,e)-pc;   

# functions given (k,b,ε)
fn_X(k,b,e) =  fn_Pi(k, e) + (1-delta) * k - b;
fn_D(next_k, next_b, k, b, e) =  fn_X(k, b, e) - next_k + q*next_b;

# Define Q(n,m,n) matrix
n = k_size * b_size * e_size  # all possible states
m = k_size * b_size           # all possible actions

# total number of possible states
s_vals = gridmake(k_grid, b_grid, e_chain.state_values)  # combination of each by value 
s_i_vals = gridmake(1:k_size, 1:b_size, 1:e_size)        # combination of each by index

a_vals = gridmake(k_grid, b_grid)
a_i_vals = gridmake(1:k_size, 1:b_size)

#################################### Standard approach
Q = zeros(n, m, n);
for next_s_i in 1:n
    for a_i in 1:m
        for s_i in 1:n

        # current state    
        e = s_i_vals[s_i, 3]   # current k, and b, does not matter for transition probabilities - enough to save e,

        # actions
        b = a_i_vals[a_i, 2]   
        k = a_i_vals[a_i, 1]   

        # future states
        next_e = s_i_vals[next_s_i, 3]
        next_b = s_i_vals[next_s_i, 2]
        next_k = s_i_vals[next_s_i, 1]

        if next_b == b && next_k == k
            Q[s_i, a_i, next_s_i] = e_chain.p[e, next_e] 
        end

        end
    end
end

R = fill(-Inf,  n, m);
for new_a_i in 1:m   
     
    # actions
    next_b = a_vals[new_a_i,2]
    next_k = a_vals[new_a_i,1]

    for s_i in 1:n

        e = s_vals[s_i, 3]    
        b = s_vals[s_i, 2]     
        k = s_vals[s_i, 1]

        d = fn_D(next_k, next_b, k, b, e)
        if d > 0
            R[s_i, new_a_i] = d  
        end

    end
end

Rhelp = [[fill("NaN", 2,3); s_vals] [transpose(a_vals) ; R]]

ddp = DiscreteDP(R, Q, discount)
@elapsed results = solve(ddp, PFI)

# Analyze results

v = results.v            # values for each state
sigma = results.sigma    # optimal policy     

reshelp = [s_vals v sigma]


##################### THE SAME PROBLEM BUT WHT STATE - ACTION PAIRS  

# R_mat - returns for each state-action pair (same as R begore)
R_mat = fill(-Inf, n, m);
@elapsed for new_a_i in 1:m  # foreach action
     
    # actions
    next_b = a_vals[new_a_i,2]
    next_k = a_vals[new_a_i,1]

    for s_i in 1:n  # foreach state

        e = s_vals[s_i, 3]    
        b = s_vals[s_i, 2]     
        k = s_vals[s_i, 1]

        d = fn_D(next_k, next_b, k, b, e)
        if d > 0
            R_mat[s_i, new_a_i] = d  
        end

    end
end

# Extract indicies for R_mat > 0
inds = Tuple.(findall(R_mat .> 0))
s_indices = first.(inds)
a_indices = last.(inds)

# Extract R_mat > 0 and arrange in a vector
mask = R_mat .> 0
R = R_mat[mask]
 
L = count(x -> x == true, R_mat .> 0)
# Construct the Q matrix
Q = spzeros(L, n);
@elapsed for next_s_i in 1:n

    # future states
    next_e = s_i_vals[next_s_i, 3]
    next_b = s_i_vals[next_s_i, 2]
    next_k = s_i_vals[next_s_i, 1]

    for sa_i in 1:L

        s_ind = s_indices[sa_i]
        a_ind = a_indices[sa_i]
            
        # current state    
        e = s_i_vals[s_ind, 3]

        b = a_i_vals[a_ind, 2]
        k = a_i_vals[a_ind, 1]

        if next_b == b && next_k == k
            Q[sa_i, next_s_i] = e_chain.p[e, next_e] 
        end
        
    end
end

@elapsed ddp = DiscreteDP(R, Q, discount, s_indices, a_indices);
@elapsed results = solve(ddp, PFI)


# interpreting results
iterations = results.num_iter
values = results.v
policy = results.sigma



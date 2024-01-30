########################################################################
######### VER 3 - ABL Qs - cash on hand - AR(1) productivity ###########
########################################################################

using LinearAlgebra, Statistics,  LaTeXStrings, Plots, QuantEcon, SparseArrays, ProgressMeter

###  State variables (x, ε)
###  Action variables (k', b')
###  ε is stochastic and engonenous 
###  x is semi-engonenous 

# grids size
x_size = 100;
k_size = 100;
b_size = 100;

# idiosyncratic shocks
e_chain = MarkovChain([0.9 0.1; 0.1 0.9], [8; 10])
e_size = length(e_chain.state_values)

# paramaters 
wage = 2.5;                # wage (exo. in the PE model)
alpha = 1/3*0.88;          # capital share
nu = 2/3*0.88;             # labour share
pc = 100;                  # participation cost
beta = 0.97;               # discount rate
delta = 0.02;              # depreciation
pdef = 0.04;               # default probability
discount = beta*(1-pdef);  # discount parameter
phi_a = 0.5;               # re-sale value of capital
phi_c = 0.5;               # variable cost of reorganization                
xi = 2*10^5;               # fixed cost of reorganization


# functions: contenporaneous optimization
fn_N(k,e) =  (nu*e*k^alpha / wage)^(1/(1-nu));      
fn_Y(k, e) =  e*k^alpha*fn_N(k,e)^nu;                 
fn_Pi(k, e) = (1-nu)*fn_Y(k,e)-pc;   

# functions: cash on hand and future values
fn_X(k,b,e) =  fn_Pi(k, e) + (1-delta) * k - b;
fn_Q(next_k, next_b) = (beta*((1-pdef)*next_b + pdef * min(next_b, phi_a*(1-delta)*next_k)))/next_b
fn_D(next_k, next_b, x, q) =  x - next_k + q*next_b; # here, defined through x

# liquidation decision

# Creating grids
# LIN GRIDS
k_grid = 2*range(3.8*10^4, 10^5, k_size)
b_grid = 1 .* range(-(10^5), 10^5, b_size) 

x_low = max(fn_X(k_grid[1],b_grid[end], e_chain.state_values[1]),0)
x_high = fn_X(k_grid[end], b_grid[1], e_chain.state_values[end])
x_grid = range(x_low, x_high, x_size)

#= LOG GRIDS
# k-grid
  k_grid = 10 .^ range(0, 6, k_size)
# this is a log grid between -5*10e5 and 5*10e5, but it is centered around 0
b_grid = 5 .* [sort(-10 .^ range(0, 5, div(b_size, 2))[2:end], rev = false); 10 .^ range(0, 5, div(b_size, 2) +1 )] 
=#

# Define Q(n,m,n) matrix
n = x_size * e_size  # all possible states
m = k_size * b_size  # all possible actions

# total number of possible states
s_vals = gridmake(x_grid, e_chain.state_values)  # combination of each by value 
s_i_vals = gridmake(1:x_size, 1:e_size)        # combination of each by index

a_vals = gridmake(k_grid, b_grid)
a_i_vals = gridmake(1:k_size, 1:b_size)

#################################### Standard approach
Q = zeros(n, m, n); 
@elapsed for next_e in 1:e_size
for a_i in 1:m
   for s_i in 1:n

       # productivities (indicies)
       e = s_i_vals[s_i, 2]  # enough to save e, current x does not matter, since there are no financial frictions
      
       # actions (values)
       b = a_vals[a_i, 2]   
       k = a_vals[a_i, 1]   
       
       # next period cash on hand x'(b',k',e')
       x_next = fn_X(k,b,e_chain.state_values[next_e])
       # where x falls on the grid  - interporalte on the grid
       x_index = argmin(abs.(x_next .- x_grid))

       # search for where it is on the s_i grid
       target_entry = (x_index, next_e) 
       xe_index = findall(r -> (r[1], r[2]) == target_entry, eachrow(s_i_vals))[1]

       # filling the transition matrix
       # probability of transition from e_i to next_e_j 
       p_trans = e_chain.p[e, next_e] 
       Q[s_i, a_i, xe_index] = p_trans
       
       end
   end
end

R = fill(-Inf,  n, m);
@elapsed for new_a_i in 1:m   
    
   # actions
   next_b = a_vals[new_a_i,2]
   next_k = a_vals[new_a_i,1]
   q = fn_Q(next_k, next_b)

   for s_i in 1:n

       e = s_vals[s_i, 2]     
       x = s_vals[s_i, 1]
       
       d = fn_D(next_k, next_b, x, q)

       if d >= 0
           R[s_i, new_a_i] = d  
       end

   end
end

ddp = DiscreteDP(R, Q, discount);
@elapsed results = solve(ddp, PFI)

# interpreting results
iterations = results.num_iter
values = results.v
policies = results.sigma   # optimal policy     

# summarising results
sumres = zeros(n, 11)
for i in 1:n

   pol = policies[i]
   k = a_vals[pol,1]
   b = a_vals[pol,2]
   e = s_vals[i, 2]
   x = fn_X(k, b, e)
   q = fn_Q(k,b)

   value = results.v[i]
   
   sumres[i, :] .= [s_vals[i, 1], s_vals[i, 2], k, b, fn_N(k, e), fn_Y(k, e), fn_Pi(k, e), x, fn_D(k, b, x, q), q, value]
  
end

sumres1 = sumres[1:Int(n/2), :]
sumres2 = sumres[Int(n/2)+1:end, :]

plot1 = plot(sumres1[:,1], [sumres1[:,3] sumres1[:,4] sumres1[:,8]], lw = 2,
   label = ["k_prime" "b_prime"  "x_prime"])
plot2 = plot(sumres1[:,1], [sumres1[:,5] sumres1[:,6] sumres1[:,7] sumres1[:,9]], lw = 2,
   label = ["n" "y" "Pi"  "d"])

plot3 = plot(sumres2[:,1], [sumres2[:,3] sumres2[:,4] sumres2[:,8]], lw = 2,
   label = ["k_prime" "b_prime"  "x_prime"])
plot4 = plot(sumres2[:,1], [sumres2[:,5] sumres2[:,6] sumres2[:,7] sumres2[:,9]], lw = 2,
   label = ["n" "y" "Pi"  "d"])
 
plot5 = plot(sumres1[:,1],  sumres1[:,10], lw = 2,
   label = "q")
plot6 = plot(sumres2[:,1], sumres2[:,10], lw = 2,
   label = "q")

   
plot(plot1, plot3, plot2, plot4, plot5, plot6, layout = (3, 2), size = (1000, 800))

  
##################### THE SAME PROBLEM BUT WITH STATE - ACTION PAIRS  

#=
R_mat = fill(-Inf,  n, m);
for new_a_i in 1:m  # same R re-labeled 
    
   # actions
   next_b = a_vals[new_a_i,2]
   next_k = a_vals[new_a_i,1]

   for s_i in 1:n

       e = s_vals[s_i, 2]     
       x = s_vals[s_i, 1]

       d = fn_D(next_k, next_b, x)

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
mask = R_mat .> 0;
R = R_mat[mask];

L = count(x -> x == true, R_mat .> 0)

# Construct the Q matrix

Q = spzeros(L, n);
@showprogress for sa_i in 1:L

   s_ind = s_indices[sa_i]
   a_ind = a_indices[sa_i]
   
   # current states
   e_i = s_i_vals[s_ind, 2]

   # current actions
   b = a_vals[a_ind, 2]
   k = a_vals[a_ind, 1]

   #probability of transition
   for next_e_i in 1:e_size

       # next period cash on hand x'(b',k',e')
       next_e = e_chain.state_values[next_e_i]
       x_next = fn_X(k,b,next_e)
       
       # where x falls on the grid  - interporalte on the grid
       x_index = argmin(abs.(x_next .- x_grid))

       # search for where it is on the s_i grid
       target_entry = (x_index, next_e_i) 
       xe_index = findall(r -> (r[1], r[2]) == target_entry, eachrow(s_i_vals))[1]

       # filling the transition matrix
       p_trans = e_chain.p[e_i, next_e_i] 
       Q[sa_i, xe_index] = p_trans
       
   end
   
end

@elapsed ddp = DiscreteDP(R, Q, discount, s_indices, a_indices);
@elapsed results = solve(ddp, PFI)


# interpreting results
iterations = results.num_iter
values = results.v
policies = results.sigma

=#
########################################################################
######### VER 3 - ABL Qs - cash on hand - AR(1) productivity ###########
########################################################################

using LinearAlgebra, Statistics, LaTeXStrings, Plots, QuantEcon, SparseArrays, ProgressMeter

###  State variables (x, ε)
###  Action variables (k', b')
###  ε is stochastic and engonenous 
###  x is semi-engonenous 

# grids size
x_size = 70;
e_size = 5;
k_size = 50;
b_size = 50;

# AR(1) produictivity
rho_e = 0.969;
sigma_e = 0.146;
nul_e = 1; 

e_chain = tauchen(e_size, rho_e, sigma_e, (1-rho_e)*nul_e, 2)
e_vals = e_chain.state_values
e_vals = exp.(e_vals) .+ 1

# paramaters 
wage = 1.85;               # wage (exo. in the PE model)
alpha = 1/3*0.88;          # capital share
nu = 2/3*0.88;             # labour share
pc = 100;                  # participation cost
beta = 0.97;               # discount rate
delta = 0.02;              # depreciation
pdef = 0.02;               # default probability
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
fn_D(next_k, next_b, x, q) =  x - next_k + q * next_b; # here, defined through x

# liquidation decision
fn_chi(k,val) = Int(phi_a*(1-delta)*k >= phi_c * val - xi)

# Creating grids
# LIN GRIDS - policy functions look much better if grids contain a 0 value
k_grid = range(0, 10^5, k_size)
b_grid =  [range(-10^5, 0, div(b_size, 2))[1:end-1]; range(0, 10^5, div(b_size, 2) +1 )]

x_low = fn_X(k_grid[1],b_grid[end], e_vals[1])
x_high = fn_X(k_grid[end], b_grid[1], e_vals[end])
x_grid =  [range(x_low, 0, div(x_size, 2))[1:end-1]; range(0, x_high, div(x_size, 2) +1 )]


#= LOG GRIDS
   # k-grid
   k_grid = 10 .^ range(0, 5, k_size)
   # this is a log grid between -5*10e5 and 5*10e5, but it is centered around 0
   b_grid = [0.5 .* sort(-10 .^ range(0, 5, div(b_size, 2))[2:end], rev = false); 10 .^ range(0, 5, div(b_size, 2) +1 )]

   x_low = fn_X(k_grid[1],b_grid[end], e_vals[1])
   x_high = fn_X(k_grid[end], b_grid[1], e_vals[end])
   x_grid = [0.5 .* sort(-10 .^ range(0, log10(-x_low), div(x_size, 2))[2:end], rev = false); 10 .^ range(0, log10(x_high), div(x_size, 2) +1 )] 
=# 

# Define Q(n,m,n) matrix
n = x_size * e_size  # all possible states
m = k_size * b_size  # all possible actions

# total number of possible states
s_vals = gridmake(x_grid, e_vals)  # combination of each by value 
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
       x_next = fn_X(k,b,e_vals[next_e])
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

       x = s_vals[s_i, 1]
       d = fn_D(next_k, next_b, x, q)

       if d >= 0
           R[s_i, new_a_i] = d 
       else
           R[s_i, new_a_i] = 1.6 * d 
       end

   end
end

ddp = DiscreteDP(R, Q, discount);
@elapsed results = solve(ddp, PFI)

# interpreting results
iterations = results.num_iter;
values = results.v;
policies = results.sigma;  # optimal policy     


# summarising results
sumres = zeros(n, 11)
for s_i in 1:n

    # policies
    pol = policies[s_i]
    k = a_vals[pol,1]
    b = a_vals[pol,2]
    
    q = fn_Q(k,b)
 
    # actions
    x = s_vals[s_i, 1]
    e = s_vals[s_i, 2]
    next_x = fn_X(k, b, e)   
    
    val = results.v[s_i]
   
   sumres[s_i, :] .= [s_vals[s_i, 1], s_vals[s_i, 2], k, b, fn_N(k, e), fn_Y(k, e), fn_Pi(k, e), next_x, fn_D(k, b, x, q), q, val]
  
end

sum_1 = sumres[sumres[:, 2] .== e_vals[1] , :]
sum_med = sumres[sumres[:, 2] .== e_vals[Int(median(1:e_size))] , :]
sum_3 = sumres[sumres[:, 2] .== e_vals[end] , :]


plot1 = plot(sum_1[:,1], [sum_1[:,3] sum_1[:,4] sum_1[:,8]], lw = 2,
   label = ["k_prime" "b_prime"  "x_prime"])
plot2 = plot(sum_1[:,1], [sum_1[:,5] sum_1[:,6] sum_1[:,7] sum_1[:,9]], lw = 2,
   label = ["n" "y" "Pi"  "d"])

plot3 = plot(sum_med[:,1], [sum_med[:,3] sum_med[:,4] sum_med[:,8]], lw = 2,
   label = ["k_prime" "b_prime"  "x_prime"])
plot4 = plot(sum_med[:,1], [sum_med[:,5] sum_med[:,6] sum_med[:,7] sum_med[:,9]], lw = 2,
   label = ["n" "y" "Pi"  "d"])

plot5 = plot(sum_3[:,1], [sum_3[:,3] sum_3[:,4] sum_3[:,8]], lw = 2,
   label = ["k_prime" "b_prime"  "x_prime"])
plot6 = plot(sum_3[:,1], [sum_3[:,5] sum_3[:,6] sum_3[:,7] sum_3[:,9]], lw = 2,
   label = ["n" "y" "Pi"  "d"])
 
plot7 = plot(sum_1[:,1],  sum_1[:,10], lw = 2,
   label = "q")
plot8 = plot(sum_med[:,1], sum_med[:,10], lw = 2,
   label = "q")
plot9 = plot(sum_3[:,1], sum_3[:,10], lw = 2,
   label = "q")
   
plot(plot1, plot3, plot5, plot2, plot4, plot6, plot7, plot8, plot9, layout = (3, 3), size = (1000, 800))


# extra codes
sumres = DataFrame(sumres, :auto)
CSV.write("matrix.csv", DataFrame(sumres))


using Pkg
Pkg.add("DataFrames")
using CSV, DataFrames


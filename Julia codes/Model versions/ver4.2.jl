########################################################################
############ VER 4.2 - default decision outside the loop ##############
########################################################################

using LinearAlgebra, Statistics, LaTeXStrings, Plots, QuantEcon, SparseArrays, ProgressMeter

###  State variables (x, ε)
###  Action variables (k', b')
###  ε is stochastic and engonenous 
###  x is semi-engonenous 

# grids size
x_size = 50;
e_size = 7;
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
pc = 1500;                 # participation cost
beta = 0.97;               # discount rate
delta = 0.02;              # depreciation
pdef_exo = 0.02;           # default probability
discount = beta*(1-pdef_exo);  # discount parameter
phi_a = 0.5;               # re-sale value of capital
phi_c = 0.5;               # variable cost of reorganization                
xi = 2*10^5;               # fixed cost of reorganization

# functions: contenporaneous optimization
fn_N(k,e) =  (nu*e*k^alpha / wage)^(1/(1-nu));      
fn_Y(k, e) =  e*k^alpha*fn_N(k,e)^nu;                 
fn_Pi(k, e) = (1-nu)*fn_Y(k,e)-pc;   

# functions: cash on hand and future values
fn_X(k,b,e) =  fn_Pi(k, e) + (1-delta) * k - b;
fn_Q(next_k, next_b, pdef) = (beta*((1-pdef)*next_b + pdef * min(next_b, phi_a*(1-delta)*next_k)))/next_b
fn_D(next_k, next_b, x, q) =  x - next_k + q * next_b; # here, defined through x

# liquidation decision
fn_chi(k,val) = Int(phi_a*(1-delta)*k >= phi_c * val - xi)

# Creating grids - should contain 0 in every case 
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
@elapsed for a_i in 1:m
   for s_i in 1:n

    # productivities (indicies)
       e = s_i_vals[s_i, 2]  # enough to save e, current x does not matter
      
       # actions (values)
       b = a_vals[a_i, 2]   
       k = a_vals[a_i, 1]   
       
       for next_e_i in 1:e_size 
            # next period cash on hand x'(b',k',e')
            x_next = fn_X(k,b,e_vals[next_e_i])
            # where x falls on the grid  - interporalte on the grid
            x_index = argmin(abs.(x_next .- x_grid))

            # search for where it is on the s_i grid
            xe_index = x_index + (next_e_i-1)*x_size

            # filling the transition matrix
            # probability of transition from e_i to next_e_j 
            p_trans = e_chain.p[e, next_e_i] 
            Q[s_i, a_i, xe_index] = p_trans
       end
    
   end
end

# initital (!) endogeneous default probability for each state
pdef_endo = zeros(n,3)
################
R = fill(-Inf,  n, m);
for new_a_i in 1:m   
    
   # actions
   next_b = a_vals[new_a_i,2]
   next_k = a_vals[new_a_i,1]

   for s_i in 1:n
       
        pdef = 1-((1-pdef_exo) * (1-pdef_endo[s_i,3]))
        q = fn_Q(next_k, next_b, pdef)

        x = s_vals[s_i, 1]
        d = fn_D(next_k, next_b, x, q)

        if d >= 0
           R[s_i, new_a_i] = d 
        else
           R[s_i, new_a_i] = 10^3 * d # try setting a proibitive cost
        end
    end
end

ddp = DiscreteDP(R, Q, discount);
@elapsed results = solve(ddp, PFI)

# interpreting results
iterations = results.num_iter;
values = results.v;
policies = results.sigma;  # optimal policy     


##### Default decision given (x,e)
defpol = zeros(n,3)
for s_i in 1:n

    x = s_vals[s_i, 1]
    e = s_vals[s_i, 2]
    val = values[s_i]
    
    def = val < 0 # default condition: value is smaller than 0

    defpol[s_i, :] .= [x, e, def]

end

##### Probability of default given optimal k', b' policies that are given by x, e
pdef_endo = zeros(n,3)
for s_i in 1:n
   
    # policies and actions
    pol = policies[s_i]
    k = a_vals[pol,1]
    b = a_vals[pol,2]
    
    e_i = s_i_vals[s_i, 2]

     for next_e_i in 1:e_size

        p_trans = e_chain.p[e_i, next_e_i] 
        x_next = fn_X(k,b,e_vals[next_e_i])
        x_index = argmin(abs.(x_next .- x_grid))
 
        # search for where it is on the s_i grid
        target_entry = (x_index, next_e_i) 
        xe_index = findall(r -> (r[1], r[2]) == target_entry, eachrow(s_i_vals))[1]
       
        # find the associated value given the xe_index
        def = defpol[xe_index,3]
        
        # update endogeneous default probability
        pdef_endo[s_i] = pdef_endo[s_i] + p_trans*def
 
   end

   pdef_endo[s_i, :] = [s_vals[s_i, 1], s_vals[s_i, 2], pdef_endo[s_i]]

end


# summarising results
nvar_sumres = 7
sumres = zeros(n, nvar_sumres)
for s_i in 1:n

    # states
    x = s_vals[s_i, 1]
    e = s_vals[s_i, 2]

    # policies
    pol = policies[s_i]
    k = a_vals[pol,1]
    b = a_vals[pol,2]
    def = defpol[s_i, 3]
    
    # outcomes
    pdef = 1-((1-pdef_exo) * (1-pdef_endo[s_i,3]))
    q = fn_Q(k, b, pdef)

    val = values[s_i]

    if def == 0
        sumres[s_i, :] .= [x, e, k, b, pdef, q, val]
    else
        sumres[s_i, 1] = x
        sumres[s_i, 2] = e
        sumres[s_i, 3:end] .= zeros(nvar_sumres-2)
    end


end

sum_1 = sumres[sumres[:, 2] .== e_vals[1] , :]
sum_med = sumres[sumres[:, 2] .== e_vals[Int(median(1:e_size))] , :]
sum_3 = sumres[sumres[:, 2] .== e_vals[end] , :]


plot1 = plot(sum_1[:,1], [sum_1[:,3] sum_1[:,4]], lw = 2,
   label = ["k_prime" "b_prime"])
plot2 = plot(sum_1[:,1], [sum_1[:,5] sum_1[:,6] ], lw = 2,
   label = ["pdef" "q"])

plot3 = plot(sum_med[:,1], [sum_med[:,3] sum_med[:,4]], lw = 2,
   label = ["k_prime" "b_prime"])
plot4 = plot(sum_med[:,1], [sum_med[:,5] sum_med[:,6] ], lw = 2,
   label = ["pdef" "q"])

plot5 = plot(sum_3[:,1], [sum_3[:,3] sum_3[:,4]], lw = 2,
   label = ["k_prime" "b_prime"])
plot6 = plot(sum_3[:,1], [sum_3[:,5] sum_3[:,6] ], lw = 2,
   label = ["pdef" "q"])


plot(plot1, plot3, plot5, plot2, plot4, plot6, layout = (2, 3), size = (1400, 800))




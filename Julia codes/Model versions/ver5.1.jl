########################################################################
###### VER 5.1 - endogneous defaults as extra state - action pair ######
########################################################################

using LinearAlgebra, Statistics, LaTeXStrings, Plots, QuantEcon, SparseArrays, ProgressMeter, Dates, XLSX, DataFrames

###  State variables (def, x, ε)
###  Action variables (def', k', b')
###  ε is stochastic and engonenous 
###  x is semi-engonenous 

# grids sizes - x,k,b should be even numbers!!
x_size = 50;
e_size = 11;
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

function fn_Q(next_k, next_b, pdef)
   if next_b == 0
      q = 0.97  # exact number does not matter
   else
      q= (beta*((1-pdef)*next_b + pdef * min(next_b, phi_a*(1-delta)*next_k)))/next_b
   end
end

fn_D(next_k, next_b, x, q) =  x - next_k + q * next_b; # here, defined through x

# liquidation decision
fn_chi(k,val) = Int(phi_a*(1-delta)*k >= phi_c * val - xi)

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

# total number of possible states
s_vals = [gridmake(x_grid, e_vals) zeros(n)]            # combination of each by value 
s_vals = [s_vals; [0 0 1]]

s_i_vals = [gridmake(1:x_size, 1:e_size) zeros(n)] 
s_i_vals = Int.([s_i_vals; [div(x_size,2) 1 1]])  # here, I set productivity to the lowest possible value, prly will not matter

a_vals = [gridmake(k_grid, b_grid) zeros(m)]
a_vals = [a_vals; [0 0 1]]

a_i_vals = [gridmake(1:k_size, 1:b_size) zeros(m)]
a_i_vals = Int.([a_i_vals; [1 div(b_size,2) 1]])

# adjusting the gridsize
n = n+1
m = m+1
#################################### Standard approach
Q = zeros(n, m, n); 
for a_i in 1:m
    for  s_i in 1:n 

     # productivities (indicies)
     e = s_i_vals[s_i, 2]  # enough to save e, current x does not matter, since there are no financial frictions

     # actions (values)
     next_def = a_vals[a_i, 3] 
     def = s_vals[s_i,3]
     b = a_vals[a_i, 2]   
     k = a_vals[a_i, 1]   
                 
        for next_e_i in 1:e_size  
            if  def == 0 && next_def == 0
            
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
                        x_far = x_close -1
                    end
                    
                    x_close_weight = abs(x_next - x_grid[x_far]) / (abs(x_next - x_grid[x_close]) + abs(x_next - x_grid[x_far]))
        
                    # finding the correspoing indicies     
                    xe_close = x_close + (next_e_i-1)*x_size
                    xe_far = x_far + (next_e_i-1)*x_size
                    
                    # filling the transition matrix
                    Q[s_i, a_i, xe_close] = p_trans*x_close_weight
                    Q[s_i, a_i, xe_far] = p_trans*(1-x_close_weight)
        
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
pdef_endo = zeros(n,m)
kbexq_old = zeros(n,4)
kbexq_new = fill(1,n,4)
################ 
iter = 0
@elapsed while !isequal(kbexq_old,kbexq_new)
   iter = iter +1
   println(iter)
   kbexq_old = kbexq_new

R = fill(-Inf,  n, m);
for a_i in 1:m   
          
   # actions
   next_b = a_vals[a_i,2]
   next_k = a_vals[a_i,1]
   next_def =  a_i_vals[a_i,3]

   for s_i in 1:n
      
      pdef = 1-((1-pdef_exo) * (1-pdef_endo[s_i, a_i]))
      q = fn_Q(next_k, next_b, pdef)

      def = s_vals[s_i, 3]

      if next_def == 1 
         #R[s_i, a_i] = x    # with this, it is exit
         R[s_i, a_i] = 0    # with this, it is default
      end

      if def == 1      
         R[s_i, a_i] = 0
      end

      if next_def == 0 && def == 0
         
         # prev. default and cash on hand and dividends
         x = s_vals[s_i, 1]
         d = fn_D(next_k, next_b, x, q)

         if d >= 0
            R[s_i, a_i] = d 
         # else
         #   R[s_i, a_i] = 10^3 * d # proibitive cost
         end
      end

    end
end

ddp = DiscreteDP(R, Q, discount);
@elapsed results = solve(ddp, PFI)

# interpreting results
values = results.v;
policies = results.sigma;  # optimal policy     

###################################################################
# summarising results
nvar_sumres = 10
sumres = zeros(n, nvar_sumres)
for s_i in 1:n

   # states
   x = s_vals[s_i, 1]
   e = s_vals[s_i, 2]

   # policies
   pol = policies[s_i]
   k = a_vals[pol,1]
   b = a_vals[pol,2]
   def = a_vals[pol,3]

   # cash on hand if productivity stays the same
   next_x = fn_X(k, b, e)   

   pdef = 1-((1-pdef_exo) * (1-pdef_endo[s_i,pol]))
   q = fn_Q(k, b, pdef)

   d = fn_D(k, b, x, q)
   val = values[s_i]

   sumres[s_i, :] .= [x, e, k, b, next_x, def, pdef, q, d, val]
end

   #= plotting
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

   =#

if iter % 5 == 0
   column_names = [:x, :e, :k, :b, :next_x, :def, :pdef, :q, :d, :val]
   sumres_df = DataFrame(sumres, column_names)
   time_string = "$(Dates.hour(now()))_$(Dates.minute(now()))_$(Dates.second(now()))"
   XLSX.writetable("$time_string.xlsx", sheetname="sheet1", sumres_df)
end 

###############################################################################
# Probability of default given optimal k', b' policies that are given by x, e
pdef_endo = zeros(n,m)
@elapsed for s_i in 1:n
     
  # state indicies
  # x_i = s_i_vals[s_i, 1]
  e_i = s_i_vals[s_i, 2]

  for a_i in 1:m
   
      # policies given (x,e)
      next_k = a_vals[a_i,1]
      next_b = a_vals[a_i,2]
      def = a_i_vals[a_i,3]

      for next_e_i in 1:e_size

         p_trans = e_chain.p[e_i, next_e_i] 
         x_next = fn_X(next_k,next_b,e_vals[next_e_i])
         x_index = argmin(abs.(x_next .- x_grid))
   
         # search for where it is on the s_i grid
         xe_index =  x_index + (next_e_i-1)*x_size
         
         # default policy given (x,e)
         pol = policies[xe_index]
         next_def = a_i_vals[pol, 3]
         
         # update endogeneous default probability
         pdef_endo[s_i, a_i] = pdef_endo[s_i, a_i] + p_trans*next_def
         
      end 

      # Does not really matter, since payoff won't depend on q() in default
      if def == 1
         pdef_endo[s_i, a_i] = 1
      end

   end
end


kbexq_new = sumres[:, [3, 4, 5, 7]]

end


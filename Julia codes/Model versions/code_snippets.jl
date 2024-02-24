###########################################################
## liquidation probability

liq_prob = zeros(n)
# calculating expected liquidation probabilities
for s_i in 1:n

   pol = policies[s_i]
   k = a_vals[pol,1]
   b = a_vals[pol,2]
   e = s_vals[s_i, 2]
   e_i = s_i_vals[s_i, 2]

  liq_prob[s_i] = 0

for next_e_i in 1:e_size

       p_trans = e_chain.p[e_i, next_e_i] 
       x_next = fn_X(k,b,e_vals[next_e])
       x_index = argmin(abs.(x_next .- x_grid))

       # search for where it is on the s_i grid
       target_entry = (x_index, next_e_i) 
       xe_index = findall(r -> (r[1], r[2]) == target_entry, eachrow(s_i_vals))[1]
      
       # find the associated value given the xe_index
       val = values[xe_index]
        
       liq_prob[s_i] = liq_prob[s_i] + p_trans*fn_chi(k,val) 

  end
end

sumres = zeros(n, 3)
for s_i in 1:n

   pol = policies[s_i]
   k = a_vals[pol,1]
   e = s_vals[s_i, 2]
  
   liq_prob[s_i]

   sumres[s_i, :] .= [k, e, liq_prob[s_i]]


end

sum_1 = sumres[sumres[:, 2] .== e_vals[1] , :]
sum_med = sumres[sumres[:, 2] .== e_vals[Int(median(1:e_size))] , :]
sum_3 = sumres[sumres[:, 2] .== e_vals[end] , :]

plot7 = plot(sum_1[:,1],  sum_1[:,3], lw = 2,
   label = "liq_prob")
plot8 = plot(sum_med[:,1], sum_med[:,3], lw = 2,
   label = "liq_prob")
plot9 = plot(sum_3[:,1], sum_3[:,3], lw = 2,
   label = "liq_prob")

plot(plot7, plot8, plot9, layout = (1, 3), size = (1000, 800))

# NOTES
#  - Here, the takeaway was that productivities are too spread out, the liquidation probabilities are either very close to 1 or 0
#  - I found it nice that liquidation probability depends mostly on productivity and not on the current cash on hand. This means that even cash poor firms will be reorganized a lot. 
#  - Not strictly related but lower dispersion of productivities might result in a nicer statitionary distribution of capital


####################################################################
##  calculating stationary distributions

dist = stationary_distributions(results.mc)[1]

# summarising the distributions across all productivities
x_dist = zeros(x_size)
for s_i in 1:e_size
    x_dist = x_dist .+ dist[1 + (s_i-1)*x_size: s_i*x_size]
end

x_dist = bar(x_dist, xlabel="Index", ylabel="Values")

bar(x_dist, xlabel="Index", ylabel="Values") |> (p -> xlims!(p, 30, 40)) |> display



using Optim

# Define your function
function my_function(x, tau)
    # Your function definition here
    return -x^2 + tau * x
end

# Set the bounds for tau
lower_bound = 0.0
upper_bound = 1.0

# Set an initial guess for tau within the bounds
initial_tau_guess = 0.5

# Define a function that takes a parameter vector [tau] and returns the negative of your function
obj_function(tau) = -my_function(2.0, tau[1])

# Use the optimize function with bounds to find the maximum
result = optimize(obj_function, [initial_tau_guess], LBFGS(), Optim.Options(g_tol=1e-6), lower=lower_bound, upper=upper_bound)

# Extract the optimal tau value
optimal_tau = result.minimizer[1]

# Evaluate the function at the optimal tau value
max_value = -result.minimum

println("Optimal tau: $optimal_tau")
println("Maximum value: $max_value")




# Optimal CFL reliance and interest rate
function fn_Tau_Q(pdef, pliq, Pi_liq, Pi_reo, next_b, tau_vec)
   q_tau = zeros(length(tau_vec))

   if next_b == 0
       fill!(q_tau, 0.97)
   else
       q_tau .= (beta ./ next_b) .* ((1 .- pdef) .* next_b .+
            pdef .* pliq .* min.(next_b, (1 .- tau_vec) .* Pi_liq .+ tau_vec .* kappa .* Pi_liq) .+
            pdef .* (1 .- pliq) .* min.(next_b, (1 .- tau_vec) .* Pi_liq .+ tau_vec .* Pi_reo))
   end

   q, tau_index = findmax(q_tau)

   if q ≈ minimum(q_tau)
       tau = 1 - pliq
   else
       tau = tau_vec[tau_index]
   end

   return q, tau
end




# Example data (replace with actual values)
pdef = 0.02

pliq = 0.25

Pi_liq = 49000
Pi_reo = 68957
next_b = 65306
tau_vec = 0.0:0.1:1.0

function plot_Tau_Q(pdef, pliq, Pi_liq, Pi_reo, next_b, tau_vec)

    q_tau = zeros(length(tau_vec))
 
    if next_b == 0
        fill!(q_tau, 0.97)
    else
        q_tau .= (beta ./ next_b) .* ((1 .- pdef) .* next_b .+
             pdef .* min.(next_b, pliq .* ((1 .- tau_vec) .* Pi_liq .+ tau_vec .* kappa .* Pi_liq) .+ 
             (1 .- pliq) .* ((1 .- tau_vec) .* Pi_liq .+ tau_vec .* Pi_reo)))
    end

    plot(tau_vec, q_tau, label="q vs tau", xlabel="tau", ylabel="q", legend=:topleft)
    title!("q vs tau") |> display

 end

plot_Tau_Q(pdef, pliq, Pi_liq, Pi_reo, next_b, tau_vec)


function fn_Tau_Q(pdef, pliq, Pi_liq, Pi_reo, next_b, tau_vec)
    q_tau = zeros(length(tau_vec))

    if next_b == 0
        fill!(q_tau, 0.97)
    else
        q_tau .= (beta ./ next_b) .* ((1 .- pdef) .* next_b .+
             pdef .* min.(next_b, pliq .* ((1 .- tau_vec) .* Pi_liq .+ tau_vec .* kappa .* Pi_liq) .+ 
             (1 .- pliq) .* ((1 .- tau_vec) .* Pi_liq .+ tau_vec .* Pi_reo)))
    end

    q, tau_index = findmax(q_tau)

    if isapprox(q, minimum(q_tau), atol=0.002) # == won't work here, there will always be a small numerical diff
        tau = Pi_reo > Pi_liq ? 1 : 0
    else
        tau = tau_vec[tau_index]
    end

    return q, tau
end
fn_Tau_Q(pdef, pliq, Pi_liq, Pi_reo, next_b, tau_vec)
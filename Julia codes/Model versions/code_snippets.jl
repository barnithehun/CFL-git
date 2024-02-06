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

########################################################################
############################### ENTRIES ################################
########################################################################

using LinearAlgebra, Statistics, LaTeXStrings, Plots, QuantEcon, Distributions

# mapping the productivity process
e_vals = e_chain.state_values
stationary_distributions(e_chain)

# entrant log-productivity distribution
#   here, set to be equal to the stationary distribution of e
#   alternatively, set it independently similar to e_chain (they you prly have to scale up the probabilities)
e_entry  = reduce(+,stationary_distributions(e_chain))

# entrant X distribution
x_vals = unique(SumPol[:, 1])
x_dist = 200*LogNormal(0.4,0.4)-300


x_entry = zeros(length(x_vals))
for i = eachindex(x_vals)
    if i == 1
        x_entry[i] = cdf(x_dist, x_vals[i])
    else
        x_entry[i] = cdf(x_dist, x_vals[i]) - cdf(x_dist, x_vals[i-1])
    end
end

# (x,e) are independent, the joint of the two distribution is their product
xe_entry = kron(x_entry, e_entry)

# map the entry probabilities to values
Ve = transpose(xe_entry) * (SumPol[1:(end-1),end]-SumPol[1:(end-1),1])


# Generate a range of x values for plotting
x_values = -400:1:1000
# Calculate the PDF values at each x value
pdf_values = pdf(x_dist, x_values)
# Plot the lognormal distribution
plot(x_values, pdf_values, label="Lognormal PDF", xlabel="x", ylabel="Probability Density", legend=:topleft)




#=
Take cdf and subtract consequtive values to get probabilities - work with indicies, then you dont need to deal with exponents
Do the same for x-es but there you need to deal with negative values
    - take a lognormal distribution and shift it with the entrycost

Then, you have two probability distributions of random variables
 they are independent so their joint is a simple multiplication 
 - I think if you map the multiplication on the xe grid it should add up to one

Map the values in the sumres to the indicies

=#



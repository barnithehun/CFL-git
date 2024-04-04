##########################################
####### CURRENTLY UNUSED RESULTS #########
##########################################

###### Exit and Entry decision visualised ######
# entry means not exiting after the productivity and x draw - so they are the same decision
ExitPol = fill(NaN, e_size, 1)
for s_i in 1:n 
    
    e_i = Int(floor( (s_i-1) / x_size) + 1)  
    x_i = (s_i % x_size == 0) ? x_size : (s_i % x_size)

    if x_i != 1 && xpol[s_i] != 1 && xpol[s_i-1] == 1
        ExitPol[e_i,1] = SumPol[x_i,1]
    end
end
plot(string.(round.(exp.(e_chain.state_values))), ExitPol, title=var, xrotation=45, legend=false, linewidth=3)

#### Firm distribution evolution - LLN ####
# how firm distribution evolves over time given some starting distribution f_dist
function xe_evolution(SumPol, Fmat, periods)

    n = size(SumPol,1)-1
    xpol = SumPol[1:n,6] + SumPol[1:n,7] # i do not account for exo exits
    xpol_mat = Matrix(I,n,n) - Diagonal(xpol)
    Mmat = Fmat*xpol_mat                 # transition matrix if entry

    f_dist = fill(1.0, n, 1)
    f_dist = (Mmat^periods)*f_dist

    # stationary x distribution
    x_dist = zeros(x_size)
    for s_i in 1:e_size
        x_dist += f_dist[1 + (s_i-1)*x_size : s_i*x_size]
    end

    # stationary e distribution
    e_dist = zeros(e_size)
    for s_i in 1:e_size
        e_dist[s_i] = sum(f_dist[1 + (s_i-1)*x_size : s_i*x_size])
    end

    xedist = plot(bar(string.(round.(exp.(e_chain.state_values))), e_dist, title = "e_dist"),
        bar(string.(round.(unique(SumPol[:, 1])./1000)),  x_dist, title = "x_dist"), 
        layout=(2,1), size=(1200, 800))

    return (xedist)

end
xe_evolution(SumPol, Fmat, 10)
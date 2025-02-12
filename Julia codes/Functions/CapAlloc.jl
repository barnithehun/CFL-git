function CapAlloc(SumPol,SumPol0,mu,mu0)
        
    x_size, e_size, _, _ = gridsize()
    n = size(SumPol,1)-1
    # this contains firms that exit immidiately which are not counted in the mode

    totK =  transpose(mu[1:n])*SumPol[1:n,3]
    totK0 =  transpose(mu0[1:n])*SumPol0[1:n,3]
    
    evec = unique(SumPol[1:n,2]) 

    CapShare = zeros(length(evec))
    CapShare0 = zeros(length(evec))

    for ind in 1:n

        e_i = div(ind-1, x_size) + 1

        CapShare[e_i] += SumPol[ind,3]*mu[ind]
        CapShare0[e_i] += SumPol0[ind,3]*mu0[ind]


    end

    CapShare =  CapShare ./ totK
    CapShare0 = CapShare0 ./ totK0

     # CapShare =  cumsum(CapShare ./  totK)
     # CapShare0 = cumsum(CapShare0 ./ totK0)

    # Logical mask: where both CapShare and CapShare0 are non-zero
    non_zero_indices = (CapShare .!= 0) .& (CapShare0 .!= 0)

    x_axis = 1:e_size

    # Filter x_axis, CapShare, and CapShare0
    filtered_x = x_axis[non_zero_indices]
    filtered_CapShare = CapShare[non_zero_indices]
    filtered_CapShare0 = CapShare0[non_zero_indices]

    # Plot the filtered data
    plott = plot(filtered_x, filtered_CapShare, label="baseline", color=:blue, linewidth=3)
    plot!(filtered_x, filtered_CapShare0, label="lowcost", color=:red, linewidth=3)


end

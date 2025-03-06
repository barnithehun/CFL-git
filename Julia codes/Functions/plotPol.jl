function plotPol(SumPolabl, SumPolcfl, e_val)
        
    SPabl = copy(SumPolabl)
    SPcfl = copy(SumPolcfl)
    xvals = string.(Int.(round.(unique(SPabl[:, 1]))))

    n = size(SPabl, 1)
    for s_i in 1:n
    
        SPabl[s_i, 6] + SPabl[s_i, 7] == 1 ? SPabl[s_i, :] .= NaN : SPabl[s_i, :]
        SPcfl[s_i, 6] + SPcfl[s_i, 7] == 1 ? SPcfl[s_i, :] .= NaN : SPcfl[s_i, :]
        
    end
    
    x_size = gridsize().x_size
    e_ind = e_val*x_size
    
    plott1 = plot(xvals, SPabl[e_ind-x_size+1:e_ind, 3], label="K", linewidth=3, alpha=0.8, xrotation=45)
    plot!(xvals, SPcfl[e_ind-x_size+1:e_ind, 3], label="K0", linewidth=3, linestyle=:dash, alpha=0.8)
    
    plott2 =plot(xvals, SPabl[e_ind-x_size+1:e_ind, 4], label="B", linewidth=3, alpha=0.8, xrotation=45)
    plot!(xvals, SPcfl[e_ind-x_size+1:e_ind, 4], label="B0", linewidth=3, linestyle=:dash, alpha=0.8)
    
    plott3 =plot(xvals, SPabl[e_ind-x_size+1:e_ind, 5], label="X", linewidth=3, alpha=0.8, xrotation=45)
    plot!(xvals, SPcfl[e_ind-x_size+1:e_ind, 5], label="X0", linewidth=3, linestyle=:dash, alpha=0.8)
    
    plott5 =plot(xvals, SPabl[e_ind-x_size+1:e_ind, 9], label="q", linewidth=3, alpha=0.8, xrotation=45)
    plot!(xvals, SPcfl[e_ind-x_size+1:e_ind, 9], label="q0", linewidth=3, linestyle=:dash, alpha=0.8)
    
    plott6 =plot(xvals, SPabl[e_ind-x_size+1:e_ind, 11], label="y", linewidth=3, alpha=0.8, xrotation=45)
    plot!(xvals, SPcfl[e_ind-x_size+1:e_ind, 11], label="y0", linewidth=3, linestyle=:dash, alpha=0.8)
    
    plott7 =plot(xvals, SPabl[e_ind-x_size+1:e_ind, 18]-SPabl[e_ind-x_size+1:e_ind, 13], label="Continuation Value", linewidth=3, alpha=0.8, xrotation=45)
    plot!(xvals, SPcfl[e_ind-x_size+1:e_ind, 18]-SPcfl[e_ind-x_size+1:e_ind, 13], label="Continuation Value0", linewidth=3, linestyle=:dash, alpha=0.8)
    
    return(plot(plott1, plott2, plott3, plott5, plott6, plott7, layout=(3,3), size=(1200, 1100)))
    
                
end


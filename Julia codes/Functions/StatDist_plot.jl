function StatDist(binnum, var::Char, SumPol)

    char_to_number = Dict('x' => 1, 'e' => 2, 'k' => 3, 'b' => 4, 'q' => 9, 'l' => 10, 'y' => 11, 'p' => 12, 'd' => 13, 't' => 17, 'v' => 18)
    varnum = get(char_to_number, var, 0)
    max = maximum(SumPol[:, varnum])

    bins = [0; exp.(range(log(10), log(max+1), binnum-1))]
    binfill = zeros(binnum, 1)
    
    n = size(SumPol,1)-1
    for s_i = 1:n

        val = SumPol[s_i,varnum] 
        val_close = argmin(abs.(val .- bins))

        if val_close != binnum
                
            if val > bins[val_close]
                binfill[val_close + 1] += mu[s_i]
            else
                binfill[val_close] += mu[s_i]
            end

        else
            binfill[val_close] += mu[s_i]
        end        
    end

    PDF = binfill ./ sum(binfill)
    CDF = cumsum(PDF, dims = 1)

    return ( PDF, CDF, bins ) 
    
end

function plotPDF(binnum, var::Char, SumPol)

    PDF, _, bins = StatDist(binnum, var::Char, SumPol)
    bin_labels = string.(round.(bins ./ 1000, digits=2))
    plot = bar(bin_labels, PDF, title=var, xrotation=45)
    return ( plot )

end

function plotCDF(binnum, var::Char, SumPol)

    _, CDF, bins = StatDist(binnum, var::Char, SumPol)
    bin_labels = string.(round.(bins ./ 1000, digits=2))
    plott = plot(bin_labels, CDF, title=var, xrotation=45, legend=false, linewidth=3)
    return ( plott )

end

# XE distribution
function plotXE(SumPol, mu, e_chain)

    x_size, e_size, _, _ = gridsize()
        x_dist = zeros(x_size)
    for s_i in 1:e_size
        x_dist += mu[1 + (s_i-1)*x_size : s_i*x_size]
    end

    # stationary e distribution
    e_dist = zeros(e_size)
    for s_i in 1:e_size
        e_dist[s_i] = sum(mu[1 + (s_i-1)*x_size : s_i*x_size])
    end

    return (plot(bar(string.(round.(exp.(e_chain.state_values))), e_dist, title = "e_dist"),
            bar(string.(round.(unique(SumPol[:, 1])./1000)),  x_dist, title = "x_dist"), 
            layout=(2,1), size=(1200, 800)))

end



# Production distribution
function ProdDist(SumPol, mu, SumPol0, mu0)

    totY =  transpose(mu)*SumPol[:,11]
    x_size, e_size, _, _ = gridsize()
    pdf_prod = zeros(e_size)
    cdf_prod = zeros(e_size)
    
    for e_ind in 1:e_size
        
        pdf_intv = (e_ind-1)*x_size+1:e_ind*x_size
        pdf_prod[e_ind] = (transpose(mu[pdf_intv]) * SumPol[pdf_intv,11]) / totY
    
        cdf_intv = 1:e_ind*x_size
        cdf_prod[e_ind] = (transpose(mu[cdf_intv]) * SumPol[cdf_intv,11]) / totY
        
    end

    totY = transpose(mu0)*SumPol0[:,11]
    x_size, e_size, _, _ = gridsize()
    pdf_prod0 = zeros(e_size)
    cdf_prod0 = zeros(e_size)

    for e_ind in 1:e_size
        
        pdf_intv = (e_ind-1)*x_size+1:e_ind*x_size
        pdf_prod0[e_ind] = (transpose(mu0[pdf_intv]) * SumPol0[pdf_intv,11]) / totY

        cdf_intv = 1:e_ind*x_size
        cdf_prod0[e_ind] = (transpose(mu0[cdf_intv]) * SumPol0[cdf_intv,11]) / totY
        
    end
    
    x_axis = 1:e_size
    plott1 = plot(x_axis, pdf_prod0, label="Production0", linewidth=3)
              plot!(x_axis, pdf_prod, label= "Production", linestyle=:dash, linewidth=3)


    plott2=  plot(x_axis, cdf_prod0, label="Cumulative Production0", linewidth=3)
              plot!(x_axis, cdf_prod, label= "Cumulative Production", linestyle=:dash, linewidth=3)
    
    return ( plot(plott1, plott2, layout=(1,2), size=(900, 400)) )     
    
end

# CFL reliance distribution
function plotTauDist(SumPol)


    bins = 0:0.1:1
    binfill = zeros(length(bins), 1)
    
    n = size(SumPol,1)-1
    for s_i = 1:n

        val = SumPol[s_i,17] 
        val_close = argmin(abs.(val .- bins))

        if SumPol[s_i, 3] != 0
            if val_close != binnum
                    
                if val > bins[val_close]
                    binfill[val_close + 1] += mu[s_i]
                else
                    binfill[val_close] += mu[s_i]
                end

            else
                binfill[val_close] += mu[s_i]
            end   
        end     
    
    end

    PDF = binfill ./ sum(binfill)
    bin_labels = string.(bins)
    plot = bar(bin_labels, PDF, title = "Histogram of CF-reliances", xrotation=45)
    return ( plot )
    
end
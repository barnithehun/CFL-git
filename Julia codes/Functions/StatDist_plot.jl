function StatDist(binnum, var::Char, SumPol)

    char_to_number = Dict('k' => 3, 'b' => 4, 'l' => 10, 'y' => 11, 'p' => 12, 'v' => 18)
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
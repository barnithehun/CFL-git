# In this version, there are productivity shocks, you could consider a setup absent of them
function dynsim(SumPol, Fmat; simn_length = 10000, e_i)

    function next_si(Fvec)
        
        r = rand()  # uniform distribution between 0-1
        cumulative_prob = 0.0
        
        for (i, prob) in enumerate(Fvec)
            cumulative_prob += prob
            if r < cumulative_prob
                return i
            end
        end 

        return length(Fvec)  # default to the last state if no transition
        
    end

    x_size, e_size, _, _ = gridsize()

    simt_length = 25
    simmat_k = fill(NaN, simn_length,simt_length);
    simmat_b = fill(NaN, simn_length,simt_length);
    simmat_x = fill(NaN, simn_length,simt_length);
    for simn in 1:simn_length 

        x_i = findall(x -> x == 0.0, unique(SumPol[:, 1]))[1]
        s_i = x_i + (e_i-1)*x_size

        for simt in 1:simt_length 

            xpol = SumPol[s_i,6] + SumPol[s_i,7] 

            if xpol != 1
                
                if simt == 1
                    simmat_x[simn,simt] = SumPol[s_i,1]
                    simmat_k[simn,simt] = SumPol[s_i,3]
                    simmat_b[simn,simt] = SumPol[s_i,4]
                else
                    # Fmat is the transposed! transition matrix   
                    Fvec = Fmat[:,s_i]
                    s_i = next_si(Fvec)
                    simmat_x[simn,simt] = SumPol[s_i,1]
                    simmat_k[simn,simt] = SumPol[s_i,3]
                    simmat_b[simn,simt] = SumPol[s_i,4]
                end

            else
                break
            end
        end
    end

    meanX = zeros(simt_length)
    meanK = zeros(simt_length) 
    meanB = zeros(simt_length)
    n_share = zeros(simt_length)
    for simt in 1:simt_length 

        meanX[simt, 1] = mean(filter(!isnan, simmat_x[:, simt]))
        meanK[simt, 1] = mean(filter(!isnan, simmat_k[:, simt]))
        meanB[simt, 1] = mean(filter(!isnan, simmat_b[:, simt]))
        n_share[simt, 1] =   (simn_length - length(filter(!isnan, simmat_x[:, simt])))/simn_length

    end

    x_axis = 1:simt_length

    plott1 = plot(x_axis, meanX, label="X", linewidth=3)
    plot!(x_axis, meanK, label="K'", linewidth=3) 
    plot!(x_axis, meanB, label="B'", linewidth=3)
    
    plott2 = plot(x_axis, n_share, label="share of exits", linewidth=3)
    
    return ( plot(plott1, plott2, layout=(1,2), size=(1200, 800))   )

end
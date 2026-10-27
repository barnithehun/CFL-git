# In this version, there are productivity shocks, you could consider a setup absent of them
function auxfun(SumPol0, Fmat0, SumPol, Fmat, e_i, var::Char, simn_length, simt_length)

    char_to_number = Dict('x' => 1, 'e' => 2, 'k' => 3, 'b' => 4, 'q' => 9, 'l' => 10, 'y' => 11, 'p' => 12, 'd' => 13, 'g' => 14, 't' => 17, 'v' => 18)

    varnum = get(char_to_number, var, 0)

    function next_si(Fvec)
        
        r = rand()  # uniform distribution between 0-1
        cumulative_prob = 0.0
        
        for (i, prob) in enumerate(Fvec)
            cumulative_prob += prob
            if r < cumulative_prob
                return i
            end
        end 

        return length(Fvec)  
        
    end

    x_size, _, _, _ = gridsize()

    simmat_0 = fill(NaN, simn_length, simt_length);
    simmat = fill(NaN, simn_length, simt_length);
    
    for simn in 1:simn_length 
    
        x_i = findall(x -> x == 0.0, unique(SumPol0[:, 1]))[1]
        s_i = x_i + (e_i-1)*x_size
    
        # Loop for 0
        for simt in 1:simt_length 
    
            xpol = SumPol0[s_i, 6] + SumPol0[s_i, 7] 
    
            if xpol != 1
                
                if simt == 1
                    simmat_0[simn, simt] = SumPol0[s_i, varnum]
                else
                    # Fmat is the transposed(!) transition matrix   
                    Fvec_0 = Fmat0[:, s_i]
                    s_i = next_si(Fvec_0)
                    simmat_0[simn, simt] = SumPol0[s_i, varnum]
                end
    
            else
                break
            end
    
        end

        x_i = findall(x -> x == 0.0, unique(SumPol[:, 1]))[1]
        s_i = x_i + (e_i-1)*x_size
    
        # Loop for cfl
        for simt in 1:simt_length 
    
            xpol = SumPol[s_i, 6] + SumPol[s_i, 7] 
    
            if xpol != 1
                
                if simt == 1
                    simmat[simn, simt] = SumPol[s_i, varnum]
                else
                    Fvec = Fmat[:, s_i]
                    s_i = next_si(Fvec)
                    simmat[simn, simt] = SumPol[s_i, varnum]
                end
    
            else
                break
            end
    
        end
    end


    meanv_0 = zeros(simt_length)
    meanv = zeros(simt_length) 
    for simt in 1:simt_length 

        meanv_0[simt, 1] = mean(filter(!isnan, simmat_0[:, simt]))
        meanv[simt, 1] = mean(filter(!isnan, simmat[:, simt]))

    end

    
    x_axis = 1:simt_length
    plott = plot(x_axis, meanv_0, label= var*"_ab", color=:blue, linewidth=3)
            plot!(x_axis, meanv, label= var*"_cf",  color=:red, linewidth=3) 


    #plott = plot(x_axis,  meanv_0 ./ meanv, label= var*"_ab", color=:blue, linewidth=3)
    #             hline!([1], linestyle=:dash, color=:gray, label="")
  
 return ( plott )

end

function dynsim2(e_i, simn_length)


    simt_length = 12
    
    p = plot(
        auxfun(SumPol0, Fmat0, SumPol, Fmat, e_i, 'k', simn_length, simt_length),
        auxfun(SumPol0, Fmat0, SumPol, Fmat, e_i, 'b', simn_length, simt_length),
        auxfun(SumPol0, Fmat0, SumPol, Fmat, e_i, 'q', simn_length, simt_length),
        auxfun(SumPol0, Fmat0, SumPol, Fmat, e_i, 'g', simn_length, simt_length),
        auxfun(SumPol0, Fmat0, SumPol, Fmat, e_i, 't', simn_length, simt_length),
        auxfun(SumPol0, Fmat0, SumPol, Fmat, e_i, 'v', simn_length, simt_length), layout=(2,3), size=(1200, 800))
    plot!(p, suptitle="Firm Lifecycle Comparision", suptitlefontsize=20)

    
end



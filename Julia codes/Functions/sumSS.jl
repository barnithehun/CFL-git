function sumSS(SumPol,Fmat,f0)
        
        mu, m, xpol = stat_dist(SumPol, Fmat, f0)
        n = size(SumPol,1)
        totalmass = sum(mu)

        # this contains firms that exit immidiately which are not counted in the model
        exitmass=transpose(mu)*xpol
        entrymass = m*(1-transpose(xpol)*f0) 
        exitshare = exitmass/totalmass 

        totK =  transpose(mu[1:n-1])*SumPol[1:n-1,3]
        totB =  transpose(mu[1:n-1])*SumPol[1:n-1,4]
        totL =  transpose(mu[1:n-1])*SumPol[1:n-1,10] # Ns = Nd
        totY =  transpose(mu[1:n-1])*SumPol[1:n-1,11]

        YtoL = totY/totL

        meanL = totL / totalmass 
        meanK = totK / totalmass 
        meanY = totY / totalmass 

        # here LIE does not work bc. Im averaging ratios - loop is more readible than the vectorized version
        avg_b2a, mu_b2a, avg_gam, mu_gam, avg_q, mu_q, avg_prod, mu_prod, avg_CFL, mu_CFL, SMEshare = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        for s_i in 1:n-1

            if SumPol[s_i,4] != 0
                avg_q += mu[s_i]/totalmass * SumPol[s_i,9]
                mu_q += mu[s_i]/totalmass
            end

            if SumPol[s_i,8] != 0
                avg_gam += mu[s_i]/totalmass * SumPol[s_i,14]
                mu_gam += mu[s_i]/totalmass
            end

            if SumPol[s_i,3] != 0 
                avg_b2a += mu[s_i]/totalmass * SumPol[s_i,4] / (SumPol[s_i, 3] + SumPol[s_i, 1])
                mu_b2a += mu[s_i]/totalmass
            end  

            if SumPol[s_i,11] != 0 && SumPol[s_i, 10] != 0
                avg_prod += mu[s_i]/totalmass * SumPol[s_i,11] ./ SumPol[s_i, 10]
                mu_prod += mu[s_i]/totalmass
            end  

            if  SumPol[s_i, 4] != 0
                avg_CFL += mu[s_i]/totalmass * SumPol[s_i,17]
                mu_CFL += mu[s_i]/totalmass
            end  

            if SumPol[s_i, 4] != 0 && (SumPol[s_i, 3]+SumPol[s_i, 1]) <= 5000 # if SumPol[s_i, 3] != 0 && (SumPol[s_i, 3]+SumPol[s_i, 1]) <= 5000 
                SMEshare += mu[s_i]/totalmass
            end
            
        end
        avg_q = avg_q / mu_q
        avg_b2a = avg_b2a / mu_b2a
        avg_gam = avg_gam / mu_gam
        avg_prod = avg_prod / mu_prod # sanity check
        avg_CFL = avg_CFL / mu_CFL
        SMEshare = SMEshare / mu_CFL # mu_CFL has the same condition in the denominator

        CFshare = (transpose(mu)*(SumPol[1:n,4] .* SumPol[1:n,17]))  /  totB

        results = zeros(15,1)
        results[:,1] = vcat(totalmass, exitmass, entrymass, exitshare, YtoL, meanL, meanK, meanY, avg_b2a, avg_q, avg_gam, avg_prod, CFshare, avg_CFL, SMEshare)
        varnames = ["totalmass", "exitmass", "entrymass", "exitshare", "YtoL", "meanL", "meanK", "meanY", "avg_b2a", "avg_q", "avg_liqprob","avg_prod", "CFshare", "avg_CFL", "SMEshare"];
        results = NamedArray(results, names=( varnames, ["values"] ) ,  dimnames=("Res", "ParamVal"))

    return ( results )

end
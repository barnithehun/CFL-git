wage = 1
phi_c = 0.4073866012405985;
zeta_Rs = 3882.0297136399195

# Cases are added manually now - if you change the number of cases, adjsut the parmater 'run'
nrun = 3

# ABL case
SumPol, e_chain, Fmat =  FirmOptim(wage, 0, zeta_Rs, zeta_Rl)
SumPolBig = fill(NaN,  size(SumPol, 1) ,  size(SumPol, 2) , nrun);
FmatBig = fill(NaN,  size(Fmat, 1) ,  size(Fmat, 2) , nrun);

SumPolBig[:,:,1] = SumPol
FmatBig[:,:,1] = Fmat

# Baseline case
SumPol, e_chain, Fmat =  FirmOptim(wage, phi_c, zeta_Rs, zeta_Rl)
SumPolBig[:,:,2] = SumPol
FmatBig[:,:,2] = Fmat

# Zeta = half case 
FirmOptim(wage, phi_c, zeta_Rs/4, zeta_Rl/4)
SumPolBig[:,:,3] = SumPol
FmatBig[:,:,3] = Fmat

# Zeta = quarter case
#SumPol, e_chain, Fmat =  FirmOptim(wage, phi_c, zeta_Rs, zeta_Rl)
#SumPolBig[:,:,4] = SumPol
#FmatBig[:,:,4] = Fmat

# Zeta = tenth case
#SumPol, e_chain, Fmat = FirmOptim(wage, phi_c, zeta_Rs, zeta_Rl)
#SumPolBig[:,:,5] = SumPol
# FmatBig[:,:,5] = Fmat


# simulates a random draw from Fvec
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



# Setting up the estimation
simn_length = 100000
simt_length = 10
e_i = 19


x_size, _, _, _ = gridsize()
x_i = findall(x -> x == 0.0, unique(SumPol[:, 1]))[1]

Bsimmat = fill(NaN, simn_length, simt_length, nrun)
Qsimmat = fill(NaN, simn_length, simt_length, nrun)
Vsimmat = fill(NaN, simn_length, simt_length, nrun)
Lsimmat = fill(NaN, simn_length, simt_length, nrun)

for run in 1:nrun  
    for simn in 1:simn_length 
        
        s_i = x_i + (e_i-1)*x_size
    
        # Loop for 0
        for simt in 1:simt_length 
    
            xpol = SumPolBig[s_i, 6, run] + SumPolBig[s_i, 7, run] 
    
            if xpol != 1

                            
                # B
                if simt == 1
                    Bsimmat[simn, simt, run] = SumPolBig[s_i, 4, run]
                else
                    # Fmat is the transposed(!) transition matrix   
                    Fvec = FmatBig[: , s_i, run]
                    s_i = next_si(Fvec)
                    Bsimmat[simn, simt, run] = SumPolBig[s_i, 4, run]
                end
                
                # Q
                if simt == 1
                    Qsimmat[simn, simt, run] = SumPolBig[s_i, 9, run]
                else
                    # Fmat is the transposed(!) transition matrix   
                    Fvec = FmatBig[: , s_i, run]
                    s_i = next_si(Fvec)
                    Qsimmat[simn, simt, run] = SumPolBig[s_i, 9, run]
                end
                

                # Val
                if simt == 1
                    Vsimmat[simn, simt, run] = SumPolBig[s_i, end, run]
                else
                    # Fmat is the transposed(!) transition matrix   
                    Fvec = FmatBig[: , s_i, run]
                    s_i = next_si(Fvec)
                    Vsimmat[simn, simt, run] = SumPolBig[s_i, end, run]
                end
                

                # Uncond. Pliq
                if simt == 1
                    Lsimmat[simn, simt, run] = SumPolBig[s_i, 14, run]
                else
                    # Fmat is the transposed(!) transition matrix   
                    Fvec = FmatBig[: , s_i, run]
                    s_i = next_si(Fvec)
                    Lsimmat[simn, simt, run] = SumPolBig[s_i, 14, run]
                end
    
            else
                break
            end
        end
    end
end

# 4 stands for the  for variables I am studying
Bmeanv = zeros(simt_length, 1, nrun) 
Qmeanv = zeros(simt_length, 1, nrun) 
Vmeanv = zeros(simt_length, 1, nrun) 
Lmeanv = zeros(simt_length, 1, nrun) 

for run in 1:nrun 
    for simt in 1:simt_length 

        Bmeanv[simt, 1, run] = mean(filter(!isnan, Bsimmat[:, simt, run]))
        Qmeanv[simt, 1, run] = mean(filter(!isnan, Qsimmat[:, simt, run]))
        Vmeanv[simt, 1, run] = mean(filter(!isnan, Vsimmat[:, simt, run]))
        Lmeanv[simt, 1, run] = mean(filter(!isnan, Lsimmat[:, simt, run]))

    end
end


x_axis = 1:simt_length
function create_plot(matrix, title_label)
    p = plot(title=title_label, xlabel="simt", ylabel="Mean Value", legend=false)
    for run in 1:nrun
        plot!(x_axis, [mean(filter(!isnan, matrix[:, t, run])) for t in x_axis], label="", lw=1, alpha=0.5)
    end
    return p
end

# Generate individual plots
p1 = create_plot(Bsimmat, "Bmeanv")
p2 = create_plot(Qsimmat, "Qmeanv")
p3 = create_plot(Vsimmat, "Vmeanv")
p4 = create_plot(Lsimmat, "Lmeanv")

# Combine plots in a grid layout
plot(p1, p2, p3, p4, layout=(2,2), size=(1200, 800))


########################## Improvements
c_e_baseline = 1079.9527380038928
zeta_Rl_vec = 0:1000:8000
zeta_Rs_vec = 0:500:4000

prodimp = zeros(length(zeta_Rl_vec), length(zeta_Rs_vec))
liqprob_SME = zeros(length(zeta_Rl_vec), length(zeta_Rs_vec))
liqprob_LE = zeros(length(zeta_Rl_vec), length(zeta_Rs_vec))
d2c_SME = zeros(length(zeta_Rl_vec), length(zeta_Rs_vec))
d2c_LE = zeros(length(zeta_Rl_vec), length(zeta_Rs_vec))
CFrel_SME = zeros(length(zeta_Rl_vec), length(zeta_Rs_vec))
CFrel_LE = zeros(length(zeta_Rl_vec), length(zeta_Rs_vec))

for (ind_l, zeta_Rl) in enumerate(zeta_Rl_vec)
    for (ind_s, zeta_Rs) in enumerate(zeta_Rl_vec)

        
        SumPol, e_chain, Fmat = FirmOptim(wage, phi_c, zeta_Rs, zeta_Rl)
        c_e, f0 = EntryValue(SumPol, e_chain) ;
        mu, m, xpol = stat_dist(SumPol, Fmat, f0);

        prodimp[ind_l, ind_s] = c_e / c_e_baseline
        liqprob_SME[ind_l, ind_s] = sumSSsme(SumPol,Fmat,f0)[10,1]
        liqprob_LE[ind_l, ind_s] = sumSSsme(SumPol,Fmat,f0)[10,2]
        d2c_SME[ind_l, ind_s] = sumSSsme(SumPol,Fmat,f0)[7,1]
        d2c_LE[ind_l, ind_s] = sumSSsme(SumPol,Fmat,f0)[7,2]
        CFrel_SME[ind_l, ind_s] = sumSSsme(SumPol,Fmat,f0)[12,1]
        CFrel_LE[ind_l, ind_s] = sumSSsme(SumPol,Fmat,f0)[12,2]

    end
end

ps = heatmap(zeta_Rs_vec, zeta_Rl_vec, (prodimp .-1) ./ 2  ,
        xlabel="zeta_Rs_vec", ylabel="zeta_Rl_vec", 
        color=:viridis, title="Productivity Improvements", size = (600,400))
savefig("productivity_heatmap.png")        

# Create the individual heatmaps
p1 = heatmap(zeta_Rs_vec, zeta_Rl_vec, liqprob_SME,
             xlabel="zeta_Rs_vec", ylabel="zeta_Rl_vec",
             color=:viridis, title="Liquidation Probability - SMEs")
p2 = heatmap(zeta_Rs_vec, zeta_Rl_vec, liqprob_LE,
             xlabel="zeta_Rs_vec", ylabel="zeta_Rl_vec",
             color=:viridis, title="Liquidation Probability - LEs")
p3 = heatmap(zeta_Rs_vec, zeta_Rl_vec, d2c_SME,
             xlabel="zeta_Rs_vec", ylabel="zeta_Rl_vec",
             color=:viridis, title="Debt-to-Collateral - SMEs")
p4 = heatmap(zeta_Rs_vec, zeta_Rl_vec, d2c_LE,
             xlabel="zeta_Rs_vec", ylabel="zeta_Rl_vec",
             color=:viridis, title="Debt-to-Collateral - LEs")
p5 = heatmap(zeta_Rs_vec, zeta_Rl_vec, CFrel_SME,
             xlabel="zeta_Rs_vec", ylabel="zeta_Rl_vec",
             color=:viridis, title="CFL Reliance - SMEs")
p6 = heatmap(zeta_Rs_vec, zeta_Rl_vec, CFrel_LE,
             xlabel="zeta_Rs_vec", ylabel="zeta_Rl_vec",
             color=:viridis, title="CFL Reliance - LEs")

plot(p1, p5, p2, p6, layout=(2,2), size = (1100,700))
savefig("combined_heatmaps.png")